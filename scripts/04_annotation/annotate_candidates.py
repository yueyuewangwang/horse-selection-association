#!/usr/bin/env python3
"""Join GWAS leads to selection regions, GFF3 genes, and optional VEP CSQ."""
from __future__ import annotations
import argparse, gzip, re
from pathlib import Path
from urllib.parse import unquote
import pandas as pd

def attrs(text):
    out={}
    for item in text.split(";"):
        if "=" in item:
            k,v=item.split("=",1); out[k]=unquote(v)
    return out
def genes(path):
    opener=gzip.open if path.suffix==".gz" else open; out={}
    with opener(path,"rt",encoding="utf-8") as h:
        for line in h:
            if line.startswith("#"): continue
            f=line.rstrip("\n").split("\t")
            if len(f)!=9 or f[2]!="gene": continue
            a=attrs(f[8]); name=a.get("Name") or a.get("gene_id","")
            if name: out.setdefault(f[0],[]).append((int(f[3]),int(f[4]),name,a.get("gene_id","")))
    return out
def main():
    p=argparse.ArgumentParser(); p.add_argument("--leads",type=Path,required=True); p.add_argument("--regions",type=Path,required=True); p.add_argument("--gff3",type=Path,required=True); p.add_argument("--out",type=Path,required=True); p.add_argument("--window",type=int,default=500_000); a=p.parse_args()
    lead=pd.read_csv(a.leads,sep="\t",dtype={"chrom":str,"CHROM":str}); reg=pd.read_csv(a.regions,sep="\t",dtype={"chrom":str,"CHROM":str}); g=genes(a.gff3)
    def col(frame,*names):
        for n in names:
            if n in frame.columns:return n
        raise ValueError(f"Missing one of {names}")
    lc,lp,lt=col(lead,"chrom","CHROM"),col(lead,"pos","POS","lead_pos"),col(lead,"trait","GWAS_Trait"); rc,rs,re=col(reg,"chrom","CHROM"),col(reg,"start","REGION_START"),col(reg,"end","REGION_END")
    snp_col = col(lead, "snp", "SNP", "lead_snp")
    rows=[]
    for x in lead.itertuples(index=False):
        chrom=str(getattr(x,lc)); pos=int(getattr(x,lp)); trait=str(getattr(x,lt)); hits=reg[(reg[rc].astype(str)==chrom)&(pd.to_numeric(reg[rs])<=pos+a.window)&(pd.to_numeric(reg[re])>=pos-a.window)]
        for y in hits.itertuples(index=False):
            start,end=int(getattr(y,rs)),int(getattr(y,re)); gh=[z for z in g.get(chrom,[]) if z[0]<=end+a.window and z[1]>=max(1,start-a.window)]
            rows.append({"trait":trait,"snp":getattr(x,snp_col),"chrom":chrom,"pos":pos,"region_start":start,"region_end":end,"genes":";".join(dict.fromkeys(z[2] for z in gh)),"gene_ids":";".join(dict.fromkeys(z[3] for z in gh if z[3]))})
    a.out.parent.mkdir(parents=True,exist_ok=True); pd.DataFrame(rows).to_csv(a.out,sep="\t",index=False)
if __name__=="__main__":main()
