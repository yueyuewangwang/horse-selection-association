#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(rMVP); library(optparse) })
opts <- parse_args(OptionParser(option_list=list(
  make_option("--manifest", type="character", help="Trait manifest"),
  make_option("--geno-root", type="character"),
  make_option("--matrix-root", type="character"),
  make_option("--ncpus", type="integer", default=1)
)))
if (is.null(opts$manifest) || is.null(opts$geno_root) || is.null(opts$matrix_root)) stop("manifest, geno-root, and matrix-root are required")
dir.create(opts$matrix_root, recursive=TRUE, showWarnings=FALSE)
manifest <- read.delim(opts$manifest, stringsAsFactors=FALSE, check.names=FALSE)
for (trait in manifest$trait) {
  input <- file.path(opts$geno_root, trait)
  output <- file.path(opts$matrix_root, trait)
  if (file.exists(paste0(output, ".geno.desc")) && file.exists(paste0(output, ".kin.desc"))) next
  MVP.Data(fileBed=input, fileKin=TRUE, filePC=TRUE, out=output, priority="speed", ncpus=opts$ncpus)
}
