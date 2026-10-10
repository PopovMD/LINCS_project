#!/usr/bin/env rscript

################ iLINCS_ALGORITHM_RUNNER.R #####################
#
# perturbagen_group  ---  name of signatures perturbagen's group (compulsory option), e.g. MEK, mTOR, PI3K...
# --project_path  ---  path to a directory with the whole project
# --out_folder  ---  NAME of out-folder, where will be located obtained results of this script. the location of this folder has already been determined by the script.
# --sng_folder --- NAME of a folder with all signatures, each signature in its own file (with strict-form table of DEGs)
#
################################################################



##########################################
# Uploading signature as a LIST OF GENES #
##########################################

library(tidyverse)
library(httr)
library(jsonlite)
library(argparse)


args <- commandArgs(trailingOnly = TRUE)

# Создаем парсер
parser <- ArgumentParser(description = "which info is to obtain")

parser$add_argument("--perturbagen_group", type = "character", default = "MEK", help = "name of signatures perturbagen's group (compulsory option), e.g. MEK, mTOR, PI3K...")
parser$add_argument("--project_path", type = "character", default = "Users/svetlana/Desktop/CMap_2025-26/", help = "path to a directory with the whole project")
parser$add_argument("--out_folder", type = "character", default = "iLINCS_results", help = "path to a directory with TDEs")
parser$add_argument("--sgn_folder", type= "character",
                    default="signatures",
                    help="The NAME of a folder with all signatures, each signature in its own file (with strict-form table of DEGs).")

# Парсим аргументы
args <- parser$parse_args()


pert_grp <- args$perturbagen_group
project_path <- args$project_path
out_folder <- args$out_folder
sgn_folder <- args$sgn_folder

sgn_directory <- paste0("/", project_path, "/", pert_grp, "_inhibitors/", sgn_folder)
results_directory <-  paste0("/", project_path, "/", pert_grp, "_inhibitors/iLINCS/", out_folder)

setwd(sgn_directory)
files <- list.files()

apiUrl="https://www.ilincs.org/api/ilincsR/findConcordancesSC"

log_file <-  paste0("/", project_path, "/", pert_grp, "_inhibitors/iLINCS/iLINCS_gss.logs")


cat(rep("\n", 5), 
    file = log_file, 
    append = T)
cat("ПОПЫТКА № 2\n\n", 
    file = log_file, 
    append = T)


for (f in files) {
  cat(paste("processing sgn:", f), 
      file = log_file, 
      append = T)
  
  setwd(sgn_directory)
  
  sgn <- read_tsv(f, col_names = F, progress = FALSE, show_col_types = FALSE)
  # Проверка минимального количества генов
  if (nrow(sgn) < 3) {
    cat(paste("СЛИШКОМ МАЛО ГЕНОВ:", nrow(sgn), "\n\n"),
        file = log_file, 
        append = T)
    next
  }
  sgn <- sgn %>% 
    select(1, 3)
  colnames(sgn) <- c('Name_GeneSymbol', 'Value_LogDiffExp')
  
  #topUpRegulatedGenes <- list(genesUp=sgn$Name_GeneSymbol[sgn$Value_LogDiffExp > 0])
  #topDownregulatedGenes <- list(genesDown=sgn$Name_GeneSymbol[sgn$Value_LogDiffExp < 0])
  genes_up <- sgn$Name_GeneSymbol[sgn$Value_LogDiffExp > 0]
  genes_down <- sgn$Name_GeneSymbol[sgn$Value_LogDiffExp < 0]
  
  body_json <- toJSON(
    list(
      mode = "UpDn",
      metadata = TRUE,
      signatureProfile = list(
        genesUp = genes_up,
        genesDown = genes_down
      )
    ),
    auto_unbox = TRUE,
    pretty = TRUE
  )
  
  try(
    req <- POST(
      "https://www.ilincs.org/api/ilincsR/findConcordancesSC",
      body = body_json,
      encode = "raw",
      content_type_json()
    )
  )
  
  cat(paste('Status_code:', req$status_code), 
      file = log_file, 
      append = T)
  #cat(paste('Status_code:', req$status_code), "\n", file = "iLINCS.log", append = TRUE)
  full_response <- content(req, "parsed")
  cat(paste("Comment:", full_response$Remark, "\n\n"), 
      file = log_file, 
      append = T)
  
  setwd(results_directory)

  try(
    concordance_table <- data.table::rbindlist(httr::content(req)$concordanceTable, use.names = TRUE, fill = TRUE)
  )
  write.table(concordance_table, 
              file = paste(sub(".txt$", "", f), "ilincs.res.tsv"), 
              sep = "\t", 
              row.names = FALSE, 
              quote = FALSE)
  
  #Задержка для соблюдения rate limits
  Sys.sleep(4)
  
}


