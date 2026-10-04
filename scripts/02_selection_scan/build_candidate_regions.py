#!/usr/bin/env python3
"""Build joint-tail selection windows, merged regions, and GFF3 genes."""
from __future__ import annotations

import argparse
import gzip
from pathlib import Path
from urllib.parse import unquote

import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--raw-root", type=Path, required=True)
    p.add_argument("--gff3", type=Path, required=True)
    p.add_argument("--out-root", type=Path, required=True)
    p.add_argument("--traits", nargs="+", required=True)
    p.add_argument("--percentile", type=float, default=5.0)
    p.add_argument("--merge-gap", type=int, default=100_000)
    p.add_argument("--gene-flank", type=int, default=100_000)
    p.add_argument("--pi-ratio", choices=("control_over_case", "case_over_control"), default="control_over_case")
    return p.parse_args()


def upper_tail(values: pd.Series, percentile: float) -> float:
    x = pd.to_numeric(values, errors="coerce").to_numpy(float)
    x = x[np.isfinite(x)]
    if not len(x):
        raise ValueError("No finite values for threshold")
    return float(np.quantile(x, 1 - percentile / 100, method="higher"))


def read_genes(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    genes = {}
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9 or f[2] != "gene":
                continue
            attrs = {}
            for item in f[8].split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    attrs[k] = unquote(v)
            name = attrs.get("Name") or attrs.get("gene_id", "")
            if name:
                genes.setdefault(f[0], []).append((int(f[3]), int(f[4]), name, attrs.get("gene_id", "")))
    return genes


def run_one(a, trait, genes):
    fst = pd.read_csv(a.raw_root / "fst" / f"{trait}.windowed.weir.fst", sep="\t", dtype={"CHROM": str})
    case = pd.read_csv(a.raw_root / "pi" / f"{trait}_case.windowed.pi", sep="\t", dtype={"CHROM": str}).rename(columns={"PI": "pi_case"})
    ctrl = pd.read_csv(a.raw_root / "pi" / f"{trait}_control.windowed.pi", sep="\t", dtype={"CHROM": str}).rename(columns={"PI": "pi_control"})
    frame = fst.merge(case[["CHROM", "BIN_START", "BIN_END", "N_VARIANTS", "pi_case"]], on=["CHROM", "BIN_START", "BIN_END", "N_VARIANTS"])
    frame = frame.merge(ctrl[["CHROM", "BIN_START", "pi_control"]], on=["CHROM", "BIN_START"])
    frame["WEIGHTED_FST"] = pd.to_numeric(frame["WEIGHTED_FST"], errors="coerce")
    num, den = (("pi_control", "pi_case") if a.pi_ratio == "control_over_case" else ("pi_case", "pi_control"))
    frame["log2_pi_ratio"] = np.log2(pd.to_numeric(frame[num], errors="coerce") / pd.to_numeric(frame[den], errors="coerce"))
    mean, sd = frame["WEIGHTED_FST"].mean(), frame["WEIGHTED_FST"].std(ddof=1)
    if not np.isfinite(sd) or sd == 0:
        raise ValueError(f"Invalid FST SD for {trait}")
    frame["ZFst"] = (frame["WEIGHTED_FST"] - mean) / sd
    frame = frame.replace([np.inf, -np.inf], np.nan).dropna(subset=["ZFst", "log2_pi_ratio"])
    zthr, pthr = upper_tail(frame.ZFst, a.percentile), upper_tail(frame.log2_pi_ratio, a.percentile)
    selected = frame[(frame.ZFst >= zthr) & (frame.log2_pi_ratio >= pthr)].sort_values(["CHROM", "BIN_START"])
    out = a.out_root / trait
    out.mkdir(parents=True, exist_ok=True)
    selected.to_csv(out / "joint_candidate_windows.tsv", sep="\t", index=False)
    merged = []
    for row in selected.itertuples(index=False):
        chrom, start, end = str(row.CHROM), int(row.BIN_START), int(row.BIN_END)
        if not merged or chrom != merged[-1]["chrom"] or start - merged[-1]["end"] > a.merge_gap:
            merged.append({"trait": trait, "chrom": chrom, "start": start, "end": end, "window_n": 1, "max_zfst": float(row.ZFst), "mean_log2_pi_ratio": float(row.log2_pi_ratio)})
        else:
            item = merged[-1]; item["end"] = max(item["end"], end); item["window_n"] += 1; item["max_zfst"] = max(item["max_zfst"], float(row.ZFst)); item["mean_log2_pi_ratio"] += float(row.log2_pi_ratio)
    for item in merged:
        item["mean_log2_pi_ratio"] /= item["window_n"]
        hits = [g for g in genes.get(item["chrom"], []) if g[0] <= item["end"] + a.gene_flank and g[1] >= max(1, item["start"] - a.gene_flank)]
        item["genes"] = ";".join(dict.fromkeys(g[2] for g in hits)); item["gene_ids"] = ";".join(dict.fromkeys(g[3] for g in hits if g[3]))
    merged_frame = pd.DataFrame(merged)
    merged_frame.to_csv(out / "merged_candidate_regions.tsv", sep="\t", index=False)
    return ({"trait": trait, "windows": len(frame), "zf_threshold": zthr, "pi_ratio_threshold": pthr,
             "joint_windows": len(selected), "merged_regions": len(merged)}, merged_frame)


def main():
    a = parse_args()
    if a.out_root.exists():
        raise FileExistsError(a.out_root)
    genes = read_genes(a.gff3)
    results = [run_one(a, trait, genes) for trait in a.traits]
    a.out_root.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([x[0] for x in results]).to_csv(a.out_root / "selection_summary.tsv", sep="\t", index=False)
    combined = pd.concat([x[1] for x in results], ignore_index=True) if results else pd.DataFrame()
    combined.to_csv(a.out_root / "candidate_regions.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
