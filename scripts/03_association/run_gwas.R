#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(rMVP); library(bigmemory); library(data.table); library(optparse) })
unlockBinding("crossprod", baseenv()); assign("crossprod", function(x,y=NULL) .Internal(crossprod(x,y)), baseenv()); lockBinding("crossprod", baseenv())
opts <- parse_args(OptionParser(option_list=list(
  make_option("--manifest", type="character"), make_option("--matrix-root", type="character"),
  make_option("--phenotype-root", type="character"), make_option("--output-root", type="character"), make_option("--ncpus", type="integer", default=1)
)))
required <- c("manifest","matrix_root","phenotype_root","output_root"); if (any(vapply(required, function(x) is.null(opts[[x]]), logical(1)))) stop("All path options are required")
manifest <- read.delim(opts$manifest, stringsAsFactors=FALSE, check.names=FALSE)
dir.create(opts$output_root, recursive=TRUE, showWarnings=FALSE)
read_aligned <- function(path, ids, keep_id=TRUE) {
  x <- read.csv(path, check.names=FALSE, stringsAsFactors=FALSE); names(x)[1] <- "Taxa"; x <- x[match(ids,x$Taxa),,drop=FALSE]
  if (anyNA(x$Taxa)) stop("Sample IDs missing from ", path)
  if (keep_id) x else { y <- x[,-1,drop=FALSE]; y[] <- lapply(y, as.numeric); y }
}
for (i in seq_len(nrow(manifest))) {
  trait <- as.character(manifest$trait[i]); prefix <- file.path(opts$matrix_root, trait); out <- file.path(opts$output_root, trait); dir.create(out, recursive=TRUE, showWarnings=FALSE)
  map <- fread(paste0(prefix,".geno.map"), header=TRUE); geno <- attach.big.matrix(paste0(prefix,".geno.desc")); kin <- attach.big.matrix(paste0(prefix,".kin.desc")); ids <- readLines(paste0(prefix,".geno.ind"))
  if (any(as.character(map$Chromosome) %in% c("X","32"))) stop("Non-autosomal marker in ", trait)
  phe <- read_aligned(file.path(opts$phenotype_root, manifest$phenotype_file[i]), ids, TRUE)
  cov <- NULL; if ("covariate_file" %in% names(manifest) && !is.na(manifest$covariate_file[i]) && nzchar(manifest$covariate_file[i])) cov <- read_aligned(file.path(opts$phenotype_root, manifest$covariate_file[i]), ids, FALSE)
  expected <- file.path(out, paste0(trait, ".MLM.csv")); if (file.exists(expected) && file.info(expected)$size > 0) next
  MVP(phe=phe, geno=geno, map=map, K=kin, nPC.MLM=10, CV.MLM=cov, ncpus=opts$ncpus, vc.method="BRENT", threshold=0.05, method="MLM", file.output=c("log","pmap","pmap.signal"), outpath=out)
}
