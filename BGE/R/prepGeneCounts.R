# Load necessary libraries
library(tidyverse)
library(tximport)
library(DESeq2)
library(sva) # For ComBat_seq
library(HTSFilter)

# Function to process data
process_expression_data <- function(file_list_path, metadata_path, expression_folder, tx2gene_path, output_path) {
  # Read file list
  file_list <- read_lines(file_list_path) # fastq read names
  
  # Read metadata
  metadata <- read.csv(metadata_path, stringsAsFactors = TRUE)
  metadata$Time <- factor(metadata$Time)
  
  # Get the name of all folders containing the quant files
  all_folders <- list.dirs(expression_folder, recursive = FALSE)
  matching_folders <- all_folders[grep("\\w+.*fam\\d", basename(all_folders))]
  files <- paste(matching_folders, "quant.sf", sep = "/")
  names(files) <- sub(".*salmon_Tcas5.2_rna_mt/(.*?)/quant\\.sf", "\\1", files)
  
  # Import transcript-level quant files and create DESeq2 object
  tx2gene <- read.csv(tx2gene_path)
  txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = FALSE, txOut = TRUE)
  ddsTC <- DESeqDataSetFromTximport(txi, colData = metadata, design = ~Condition + Time + Condition:Time)
  ddsTC <- HTSFilter(ddsTC, s.len = 60, plot = FALSE)$filteredData
  
  # Run DESeq normalization and extract counts
  ddsTC <- estimateSizeFactors(ddsTC)
  counts <- counts(ddsTC, normalized = FALSE)
  
  # Create a model matrix for the experimental design variables (without the batch effect)
  model_matrix <- model.matrix(~Condition + Time + Condition:Time, colData(ddsTC))
  
  # Define batch variable
  colData(ddsTC)$batch <- with(colData(ddsTC), paste0(Flowcell, "_", Lane))
  
  # Run ComBat for batch effect correction
  batch_corrected_counts <- ComBat_seq(counts = counts, batch = colData(ddsTC)$Flowcell, covar_mod = model_matrix)
  expressionData_raw <- t(batch_corrected_counts) # Raw for future DESeq2 differential expression analysis
  
  rlogData <- rlog(ddsTC, blind = FALSE)
  expressionData <- assay(rlogData)
  expressionData <- t(expressionData) # For WGCNA
  
  # Save the processed data
  save(expressionData, expressionData_raw, metadata, file = output_path)
}
