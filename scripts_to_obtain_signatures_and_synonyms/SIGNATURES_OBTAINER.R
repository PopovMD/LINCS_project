#!/usr/bin/env rscript

################# obtain_signatures.R #####################
#
# --tdedir  ---  path to a directory with TDEs 
# --sgndir  ---  path to a directory to use as signatures' storage
# --pval  ---  p-value range
# --lfc  ---  Log2FoldChange range
#
###########################################################


library(parallel)
library(fs)
library(argparse)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

# Создаем парсер
parser <- ArgumentParser(description = "which info is to obtain")

parser$add_argument("--pval", type = "numeric", default = 0.05, help = "Limits for p-value")
parser$add_argument("--lfc", type = "numeric", default = 1, help = "limits for Log2FoldChange value")
parser$add_argument("--tdedir", type = "character", default = "/data/popov/mTOR_inhibitors/results_TDEs", help = "path to a directory with TDEs")
parser$add_argument("--sgndir", type = "character", default = "/data/popov/mTOR_inhibitors/signatures", help = "path to a directory to use as signatures' storage")

# Парсим аргументы
args <- parser$parse_args()


pval <- args$pval
lfc <- args$lfc
TDE_dir <- args$tdedir
sgn_dir <- args$sgndir

obtain_signature <- function(tde, pval, abs_lfc) {

    DE_table <- read_tsv(tde)

    DE_table <- DE_table %>% 
        arrange(desc(abs(log10(pvalue))*sign(log2FoldChange)))

    siggenes <- DE_table %>%
        dplyr::filter(padj < pval, abs(log2FoldChange)>abs_lfc)
    
    signature <- siggenes %>%
        dplyr::select(gene_name, gene_id, log2FoldChange, pvalue, padj)
    
    return(signature)
    
}

#if (inh_grp == "MEK"){
#	TDE_dir <- paste0(pd, "/signatures_", inh_grp, "-inhibitors/results_TDEs")
#	sgn_dir <- paste0(pd, "/signatures_", inh_grp, "-inhibitors/signatures")
#} else {
#	TDE_dir <- paste0(pd, "/", inh_grp, "_inhibitors/results_TDEs")
#	sgn_dir <- paste0(pd, "/", inh_grp, "_inhibitors/signatures")
#}

setwd(TDE_dir)
files <- list.files()

tdes <- vector(mode="character", length=length(files))

for (i in 1:length(files)) {
    file = paste0(TDE_dir, "/", files[i])
    tdes[i] <- file
}

results <- mclapply(tdes, obtain_signature, pval, lfc, mc.cores = 1)

setwd(sgn_dir)

for (i in 1:length(files)) {
    out_file <- paste0(str_remove(files[i], fixed(".deseq2_table.tsv")), ".sgn.txt")
    write.table(
        x = results[i], 
        file = out_file,
        sep = "\t",
        row.names = F, 
        col.names = T,
        quote = F)
}








