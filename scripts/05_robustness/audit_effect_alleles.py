#!/usr/bin/env python3
"""Check that the declared effect allele is present in the PLINK BIM file."""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
def main():
    p=argparse.ArgumentParser(); p.add_argument("--bim",type=Path,required=True); p.add_argument("--leads",type=Path,required=True); p.add_argument("--out",type=Path,required=True); a=p.parse_args()
    bim=pd.read_csv(a.bim,sep=r"\s+",header=None,names=["CHR","SNP","CM","POS","A1","A2"],dtype=str); lead=pd.read_csv(a.leads,sep="\t",dtype=str)
    required={"SNP","effect_allele"};
    if not required.issubset(lead.columns): raise ValueError(f"lead table needs {sorted(required)}")
    out=lead.merge(bim[["SNP","A1","A2"]],on="SNP",how="left"); out["effect_allele_present"]=out.apply(lambda r: str(r.effect_allele) in {str(r.A1),str(r.A2)},axis=1)
    a.out.parent.mkdir(parents=True,exist_ok=True); out.to_csv(a.out,sep="\t",index=False)
    if not out.effect_allele_present.all(): raise SystemExit("One or more effect alleles are absent from BIM")
if __name__=="__main__": main()
