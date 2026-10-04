# Reproducible core workflow

This is the simplified public code repository for the reported horse
population-genomic analyses. It is organised by scientific module:

```text
scripts/
├── 00_snp_filter           genotype QC and preparation
├── 01_population_structure PCA, VCF2Dis distance, tree, ADMIXTURE
├── 02_selection_scan       FST, nucleotide diversity, joint regions
├── 03_association          autosomal rMVP mixed-model GWAS
├── 04_annotation           VEP and GWAS/selection annotation
└── 05_robustness           frequencies, haplotypes, sensitivity analyses
```

The workflow contains no server-specific paths, private sample metadata, large
VCFs, PLINK matrices, reference genome, or VEP cache. Set local paths in
`config/config.env`; the committed template uses only portable example paths.

## Quick start

```bash
cp config/config.example.env config/config.env
# Edit config/config.env and the input manifests.
source config/config.env
bash workflow/run_stage.sh --help
```

The scientific meaning and input contract are in [docs/ANALYSIS_MAP.md](docs/ANALYSIS_MAP.md)
and [docs/DATA_LAYOUT.md](docs/DATA_LAYOUT.md). The organisation follows the
modular public-repository convention used by projects such as
[CAU-T2T-Goat](https://github.com/Wuhui2024/CAU-T2T-Goat), where each major
analysis has its own directory and module purpose is explained in the README.

## Main stages

| Stage | Purpose | Command |
|---|---|---|
| QC | exceptional-sample removal, missingness filtering, post-PCA filtering, LD pruning | `bash workflow/run_stage.sh qc` |
| Structure | PCA, VCF2Dis distance, FastME BIONJ tree, ADMIXTURE K comparison | `bash workflow/run_stage.sh structure` |
| Selection | weighted FST, windowed pi, joint-tail regions, GFF3 genes | `bash workflow/run_stage.sh selection` |
| Association | chromosome-X removal, kinship/PC matrices, rMVP MLM, locus summary | `bash workflow/run_stage.sh association` |
| Annotation | offline VEP and physical GWAS/selection integration | `bash workflow/run_stage.sh annotation` |
| Robustness | allele frequencies, phased regions, candidate and selection sensitivity | `bash workflow/run_stage.sh robustness` |

The stages are separate because they have different inputs and computational
costs. Inspect intermediate tables before moving to the next stage. Main
wrappers do not silently overwrite an existing output directory.

## Analysis conventions

- Chromosome X is excluded before association marker counts, kinship, PCs,
  thresholds, significant markers, loci, and integration tables are created.
- Environmental covariates are supplied explicitly through the GWAS trait
  manifest; they are never inferred automatically.
- Selection thresholds are calculated within the current window universe. The
  output records window size, step, percentile, merge gap, and pi-ratio
  orientation.
- VCF2Dis calculates pairwise p-distance, FastME reconstructs the BIONJ tree,
  and ADMIXTURE Q components are descriptive model labels.
- GWAS/selection physical proximity is candidate evidence, not formal
  colocalisation or causal proof.
- Robustness scripts operate on declared candidate, bin, region, and scenario
  manifests. They do not add targeted genes or exploratory exclusions.

## Dependencies

PLINK 2, PLINK 1.9, vcftools, BCFtools, VCF2Dis, FastME, ADMIXTURE, Beagle,
and VEP with a compatible offline cache are required for the relevant stages.
R scripts use rMVP/bigmemory/data.table/optparse where indicated. Python
dependencies are listed in [requirements-python.txt](requirements-python.txt).
