#!/usr/bin/env rscript

################ iLINCS_ALGORITHM_RUNNER.R #####################
#
# perturbagen_group  ---  name of signatures perturbagen's group (compulsory option), e.g. MEK, mTOR, PI3K...
# --project_path  ---  path to a directory with the whole project
# --out_folder  ---  NAME of out-folder, where will be located obtained results of this script. the location of this folder has already been determined by the script.
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

# Парсим аргументы
args <- parser$parse_args()


pert_grp <- args$perturbagen_group
project_path <- args$project_path
out_folder <- args$out_folder


sgn_directory <- paste0("/", project_path, "/", pert_grp, "_inhibitors/signatures")
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


#ilincsUpDnConnectedSignatures <- data.table::rbindlist(httr::content(req)$concordanceTable, use.names = TRUE, fill = TRUE)

#head(ilincsUpDnConnectedSignatures)




# uncomment the next section if tou want to try to convert symbol-ids into entrez-ids 
# and try out iLINCS algorithm with entrez-ids. 
# There may be some directories missing, so be careful.

'''

##### TRYING ENTREZID #####
## смотрим как пройдет конвертация из ensembl в entrez

library(org.Hs.eg.db)
library(AnnotationDbi)

na_scores_ens <- vector(mode = "numeric", length=length(sgns))
n_of_rows <- vector(mode = "logical", length = length(sgns))
poryadok_check <- vector(mode = "logical", length = length(sgns))

for (i in 1:length(sgns)) {
  setwd(paste0("/", project_path, "/", pert_grp, "_inhibitors/signatures"))
  
  sgn <- read_tsv(sgns[i], col_names = F)
  print(paste("PROCESSING SGN №", i))
  if (nrow(sgn) != 0) {
    
    try(gene_ens <- sgn[[2]])
    
    # 3. Convert Symbols to Entrez IDs
    try(
      entrez_ens <- select(org.Hs.eg.db, 
                             keys = gene_ens, 
                             columns = c("ENTREZID", "GENENAME"), 
                             keytype = "ENSEMBL")
    )
    #na_scores_alias_1[i] <- sum(is.na(entrez_alias[["ENTREZID"]])) / length(entrez_alias[["ENTREZID"]])
    
    df <- entrez_ens %>% distinct(entrez_ens[[1]], .keep_all = TRUE)
    df <- df[,-4]
    
    # df <- df[!duplicated(df[, 1]), ]
    
    # проверка на повторения айдишников
    n_of_rows[i] <- nrow(sgn) == nrow(df)
    "
    if (n_of_rows[i] == FALSE) {
      sgn <- sgn %>% distinct(sgn[[2]], .keep_all = TRUE)
    }
    "
    df <- df %>% 
      mutate("SYMBOL" = sgn[[1]], "l2fc" = sgn[[3]], .before = 1)
    # проверка отсутствия нарушений в порядке айдишников
    poryadok_check[i] <- (sum(df[[1]] == sgn[[1]]) == nrow(sgn))
    
    na_scores_ens[i] <- sum(is.na(df[["ENTREZID"]])) / length(df[["ENTREZID"]])
    
    # сохраняем новую сигнатуру
    setwd(paste0("/", project_path, "/", pert_grp, "_inhibitors/signatures_extended(entrezid)"))
    write_tsv(df, sgns[i])
    
  } else {
    n_of_rows[i] <- TRUE
    poryadok_check[i] <- TRUE
    na_scores_ens[i] <- 0
  }
}

sum(n_of_rows) == length(n_of_rows)
n_of_rows

sum(poryadok_check) == length(poryadok_check)
poryadok_check

print(na_scores_ens)
max(na_scores_ens)
which(na_scores_ens == max(na_scores_ens))

mean(na_scores_ens) # 0.069 - без дубликатов
mean(na_scores_symbol) # 0.11
mean(na_scores_alias) # 0.05







##### POSTING ENTREZIDs #####


library(tidyverse)
library(httr)
library(jsonlite)


sgn_directory <- paste0("/", project_path, "/", pert_grp, "_inhibitors/signatures_extended(entrezid)")
results_directory <- paste0("/", project_path, "/", pert_grp, "_inhibitors/iLINCS/iLINCS_results_EntrezID")

setwd(sgn_directory)
files <- list.files()

apiUrl="https://www.ilincs.org/api/ilincsR/findConcordancesSC"

log_file <- paste0("/", project_path, "/", pert_grp, "_inhibitors/iLINCS/iLINCS_gss.logs")


cat(rep("\n", 5), 
    file = log_file, 
    append = T)
cat("ПОПЫТКА № 3 (EnsemblID -> EntrezID)\n\n", 
    file = log_file, 
    append = T)


for (f in files) {
  cat(paste("processing sgn:", f), 
      file = log_file, 
      append = T)
  
  setwd(sgn_directory)
  
  sgn <- read_tsv(f, progress = FALSE, show_col_types = FALSE)
  # Проверка минимального количества генов
  if (nrow(sgn) < 3) {
    cat(paste("СЛИШКОМ МАЛО ГЕНОВ:", nrow(sgn), "\n\n"),
        file = log_file, 
        append = T)
    next
  }
  sgn <- sgn %>% 
    dplyr::select(4, 2)
  colnames(sgn) <- c("Gene_EntrezID", "Value_LogDiffExp")
  
  #topUpRegulatedGenes <- list(genesUp=sgn$Name_GeneSymbol[sgn$Value_LogDiffExp > 0])
  #topDownregulatedGenes <- list(genesDown=sgn$Name_GeneSymbol[sgn$Value_LogDiffExp < 0])
  genes_up <- sgn$Gene_EntrezID[sgn$Value_LogDiffExp > 0]
  genes_down <- sgn$Gene_EntrezID[sgn$Value_LogDiffExp < 0]
  
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
  
  cat(paste("Status_code:", req$status_code), 
      file = log_file, 
      append = T)
  #cat(paste("Status_code:", req$status_code), "\n", file = "iLINCS.log", append = TRUE)
  full_response <- httr::content(req, as = "parsed")
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

'''


