#!/usr/bin/env python3
"""Calculate group-stratified frequencies for a declared lead-SNP table.

Inputs are generic: a PLINK prefix, a phenotype table (IID/value), a bin table
(trait/label/lower/upper), and a lead table (trait/SNP/effect_allele).
"""
from __future__ import annotations
import argparse, subprocess, tempfile
from pathlib import Path
import pandas as pd

def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument("--bfile",type=Path,required=True); p.add_argument("--phenotype",type=Path,required=True); p.add_argument("--bins",type=Path,required=True); p.add_argument("--leads",type=Path,required=True); p.add_argument("--out",type=Path,required=True); p.add_argument("--plink",default="plink"); a=p.parse_args()
    lead=pd.read_csv(a.leads,sep="\t",dtype=str); required={"trait","SNP","effect_allele"}
    if not required.issubset(lead.columns): raise ValueError(f"lead table needs {sorted(required)}")
    phe=pd.read_csv(a.phenotype,sep="\t",dtype=str); bins=pd.read_csv(a.bins,sep="\t")
    if not {"IID","value"}.issubset(phe.columns): raise ValueError("phenotype needs IID and value")
    if not {"trait","label","lower","upper"}.issubset(bins.columns): raise ValueError("bins needs trait,label,lower,upper")
    fam=pd.read_csv(str(a.bfile)+".fam",sep=r"\s+",header=None,usecols=[1],names=["IID"],dtype=str)
    snp_file=Path(tempfile.mkstemp(prefix="lead_snps_",suffix=".txt")[1]); group_file=Path(tempfile.mkstemp(prefix="groups_",suffix=".txt")[1]);
    try:
        lead.SNP.drop_duplicates().to_csv(snp_file,index=False,header=False)
        outputs=[]
        for trait, group_bins in bins.groupby("trait",sort=False):
            vals=phe[phe.IID.isin(fam.IID)].copy(); vals["value"]=pd.to_numeric(vals.value,errors="coerce"); vals=vals.dropna()
            rows=[]
            for b in group_bins.itertuples(index=False):
                subset=vals[(vals.value>=float(b.lower))&(vals.value<float(b.upper))]
                for iid in subset.IID: rows.append(f"0 {iid} {b.label}")
            group_file.write_text("\n".join(rows)+"\n")
            prefix=Path(tempfile.mkstemp(prefix="freq_",suffix="")[1]); prefix.unlink(missing_ok=True)
            cmd=[a.plink,"--bfile",str(a.bfile),"--extract",str(snp_file),"--freq","--within",str(group_file),"--allow-no-sex","--out",str(prefix)]
            subprocess.run(cmd,check=True)
            frq=pd.read_csv(str(prefix)+".frq.strat",sep=r"\s+",engine="python")
            frq.insert(0,"trait",trait); outputs.append(frq)
        a.out.parent.mkdir(parents=True,exist_ok=True); pd.concat(outputs,ignore_index=True).to_csv(a.out,sep="\t",index=False)
    finally:
        snp_file.unlink(missing_ok=True); group_file.unlink(missing_ok=True)
if __name__=="__main__": main()
