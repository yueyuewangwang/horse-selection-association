# Portable data layout

The workflow does not require a particular server layout. A local project can
use this conceptual organisation or override every location in `config.env`:

```text
project/
├── input/       genotypes, phenotypes, metadata, selection, annotation
├── reference/   genome FASTA and offline VEP cache
├── work/        intermediate files
└── results/     final tables and figures
```

Important manifests are tab-separated:

```text
# GWAS traits
trait  genotype_prefix  phenotype_file  covariate_file

# Selection groups
trait  case_file  control_file

# Candidate regions
region  chrom  start  end
```

Sample identifiers must agree across genotype, phenotype, group, PC, metadata,
and haplotype files. Large inputs and private identifiers are external.
