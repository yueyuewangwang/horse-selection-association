#!/usr/bin/env python3
"""Summarise phased VCF region sizes and sample counts."""
from __future__ import annotations
import argparse, gzip
from pathlib import Path
import pandas as pd
def main():
    p=argparse.ArgumentParser(); p.add_argument("--manifest",type=Path,required=True); p.add_argument("--vcf-root",type=Path,required=True); p.add_argument("--out",type=Path,required=True); a=p.parse_args(); rows=[]
    for r in pd.read_csv(a.manifest,sep="\t").itertuples(index=False):
        path=a.vcf_root/f"{str(r.region)}.phased.vcf.gz"; nvar=0; nsample=0
        with gzip.open(path,"rt") as h:
            for line in h:
                if line.startswith("#CHROM"): nsample=max(0,len(line.rstrip().split("\t"))-9)
                elif not line.startswith("#"): nvar+=1
        rows.append({"region":r.region,"variant_n":nvar,"sample_n":nsample})
    a.out.parent.mkdir(parents=True,exist_ok=True); pd.DataFrame(rows).to_csv(a.out,sep="\t",index=False)
if __name__=="__main__": main()
