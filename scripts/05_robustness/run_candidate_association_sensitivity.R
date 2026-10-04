#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(optparse) })
opts <- parse_args(OptionParser(option_list=list(
  make_option("--candidate-dosage", type="character"), make_option("--phenotype-manifest", type="character"),
  make_option("--pc-file", type="character"), make_option("--metadata", type="character"), make_option("--out", type="character")
)))
required <- c("candidate_dosage","phenotype_manifest","pc_file","metadata","out")
if (any(vapply(required, function(x) is.null(opts[[x]]), logical(1)))) stop("All input options are required")
dosage <- read.delim(opts$candidate_dosage, check.names=FALSE, stringsAsFactors=FALSE)
manifest <- read.delim(opts$phenotype_manifest, stringsAsFactors=FALSE, check.names=FALSE)
pcs <- read.delim(opts$pc_file, stringsAsFactors=FALSE, check.names=FALSE)
meta <- read.delim(opts$metadata, stringsAsFactors=FALSE, check.names=FALSE)
if (!"IID" %in% names(dosage) || !"IID" %in% names(pcs) || !"IID" %in% names(meta)) stop("All tables need IID")
pc_names <- intersect(paste0("PC",1:10), names(pcs)); dosage <- merge(dosage, pcs[,c("IID",pc_names)], by="IID"); dosage <- merge(dosage, meta, by="IID")
results <- list()
for (i in seq_len(nrow(manifest))) {
  trait <- manifest$trait[i]; phenotype <- read.delim(manifest$phenotype_file[i], stringsAsFactors=FALSE, check.names=FALSE)
  names(phenotype)[1:2] <- c("IID","value"); dat <- merge(dosage, phenotype[,c("IID","value")], by="IID"); dat$value <- as.numeric(dat$value)
  if ("covariate_file" %in% names(manifest) && !is.na(manifest$covariate_file[i]) && nzchar(manifest$covariate_file[i])) { cv <- read.delim(manifest$covariate_file[i], stringsAsFactors=FALSE); names(cv)[1:2] <- c("IID","covariate"); dat <- merge(dat,cv[,c("IID","covariate")],by="IID"); dat$covariate <- as.numeric(dat$covariate) }
  snps <- setdiff(names(dosage), c("IID",pc_names,names(meta),"value","covariate"))
  for (snp in snps) for (pc_set in c(2L,10L)) {
    covars <- paste(c(paste0("PC",seq_len(min(pc_set,length(pc_names)))), if("covariate" %in% names(dat)) "covariate"), collapse="+")
    formula <- as.formula(paste("value ~", snp, if(nzchar(covars)) paste("+",covars) else "")); fit <- lm(formula, data=dat); co <- summary(fit)$coefficients
    if (snp %in% rownames(co)) results[[length(results)+1L]] <- data.frame(trait=trait,SNP=snp,model=paste0("PC",pc_set),beta=co[snp,1],SE=co[snp,2],p_value=co[snp,4],sample_n=nrow(dat))
  }
}
dir.create(dirname(opts$out),recursive=TRUE,showWarnings=FALSE); write.table(do.call(rbind,results),opts$out,sep="\t",row.names=FALSE,quote=FALSE)
