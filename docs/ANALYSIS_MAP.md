# Analysis map

## `00_snp_filter`

Remove explicitly documented exceptional samples, apply individual and marker
missingness filters, apply the post-PCA sample exclusion list, and produce a
MAF-filtered LD-pruned set. This defines the analysis cohort and is not a
biological test.

## `01_population_structure`

PCA describes major genetic axes. VCF2Dis calculates the pairwise p-distance
matrix, FastME reconstructs a BIONJ tree, and ADMIXTURE estimates model-based
ancestry proportions over a declared K range with optional cross-validation.
ADMIXTURE components are descriptive and are not named source populations.

## `02_selection_scan`

For each case/control contrast, calculate weighted Weir–Cockerham FST and
windowed nucleotide diversity. Standardise FST to ZFst, intersect upper tails
of ZFst and the declared pi-ratio, merge nearby windows, and annotate regions
with GFF3. Parameters are explicit command-line inputs.

## `03_association`

Remove chromosome X before genotype conversion, rebuild autosomal kinship and
PC inputs, and fit the rMVP mixed linear model with PC1–PC10. Prespecified
environmental covariates are supplied in the trait manifest. Marker-wise
thresholds use the current autosomal marker count.

## `04_annotation`

Run offline VEP on a declared variant set and join GWAS leads to selection
regions by physical overlap/proximity. GFF3 genes and VEP consequences provide
annotation, not formal colocalisation or causal validation.

## `05_robustness`

Summarise declared alleles across phenotype/environment bins, audit allele
orientation, phase declared candidate regions, and test robustness with
prespecified candidate-association and selection scenarios. These analyses are
follow-ups and do not replace genome-wide discovery tests.
