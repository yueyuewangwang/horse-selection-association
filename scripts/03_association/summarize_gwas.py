#!/usr/bin/env python3
"""Create autosomal significant-marker and locus tables from MLM outputs."""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

def main():
    p=argparse.ArgumentParser(); p.add_argument("--manifest",type=Path,required=True); p.add_argument("--matrix-root",type=Path,required=True); p.add_argument("--result-root",type=Path,required=True); p.add_argument("--out-root",type=Path,required=True); p.add_argument("--locus-gap",type=int,default=1_000_000); a=p.parse_args()
    rows=[]; loci=[]
    manifest=pd.read_csv(a.manifest,sep="\t",dtype=str)
    for trait in manifest.trait:
        path=a.result_root/trait/f"{trait}.MLM.csv"
        if not path.is_file(): raise FileNotFoundError(path)
        n=sum(1 for _ in open(a.matrix_root/trait/f"{trait}.geno.map"))-1; threshold=0.05/n
        frame=pd.read_csv(path); pcol=next((x for x in frame.columns if x.endswith(".MLM") or x==trait),None)
        if pcol is None: raise ValueError(f"No MLM p-value column in {path}")
        frame=frame.rename(columns={pcol:"P"}); frame["P"]=pd.to_numeric(frame.P,errors="coerce"); frame["CHROM"]=frame["CHROM"].astype(str)
        if frame.CHROM.isin(["X","32"]).any(): raise ValueError(f"X chromosome present in {path}")
        sig=frame[frame.P<=threshold].copy(); sig.insert(0,"trait",trait); rows.append(sig)
        for chrom, part in sig.sort_values("POS").groupby("CHROM",sort=False):
            part=part.sort_values("POS"); start=0
            for j in range(1,len(part)+1):
                if j==len(part) or int(part.iloc[j].POS)-int(part.iloc[j-1].POS)>a.locus_gap:
                    block=part.iloc[start:j]; lead=block.loc[block.P.idxmin()]
                    loci.append({"trait":trait,"chrom":str(chrom),"start":int(block.POS.min()),"end":int(block.POS.max()),"snp_n":len(block),"lead_snp":str(lead.SNP),"lead_pos":int(lead.POS),"lead_p":float(lead.P),"lead_effect":float(lead.get("Effect",np.nan))}); start=j
    a.out_root.mkdir(parents=True,exist_ok=True); pd.concat(rows,ignore_index=True).to_csv(a.out_root/"significant_snps.tsv",sep="\t",index=False); pd.DataFrame(loci).to_csv(a.out_root/"lead_loci.tsv",sep="\t",index=False)

if __name__=="__main__": main()
