######################### Packages loading, setting up #########################

library(limma)
library(edgeR)
library(DESeq2)
library(dplyr)
library(reshape2)
library(tibble)
library(tidyverse)
library(Rtsne)
library(umap)

library(ggplot2)
library(gplots)
library(ggrepel)
library(ggplotify)
library(pheatmap)
library(heatmaply)
library(RColorBrewer)
library(grDevices)

setwd("~/Downloads/OE33_project")

########################### Reading & Preparing the data #######################

# GENES
counts_genes <- read.delim("./2023_NovoGene_RNAseq_Report/OE33_salmon.merged.gene_counts.tsv")
# TRANSCRIPTS:
counts_transcripts <- read.delim("./2023_NovoGene_RNAseq_Report/OE33_salmon.merged.transcript_counts.tsv")

dim(counts_genes) # 62703
dim(counts_transcripts) # 252301
head(counts_genes)
head(counts_transcripts)

# GENES: rownames are gene ids
rownames(counts_genes) <- counts_genes$gene_id
# TRANSCRIPTS: rownames are transcript ids
rownames(counts_transcripts) <- counts_transcripts$tx

# Removing duplicates
# GENES:
duplicated_rows <- duplicated(counts_genes$gene_name)
sum(duplicated_rows)
counts_genes <- counts_genes[!duplicated_rows, ]
# TRANSCRIPTS:
duplicated_rows <- duplicated(counts_transcripts$tx)
sum(duplicated_rows)
counts_transcripts <- counts_transcripts[!duplicated_rows, ]

############################## Creating a DESeq2 dataset #######################

# Creating a design matrix - make sure to CHECK THE ORDER of our samples!
condition <- c("Control_AB1", "Control_AB1", "Control_AB1",
               "Control_AE21", "Control_AE21", "Control_AE21",
               "METTL1_KO_AB16", "METTL1_KO_AB16", "METTL1_KO_AB16",
               "METTL1_KO_AE4", "METTL1_KO_AE4", "METTL1_KO_AE4")

# Saving names
names_genes <- counts_genes[,1:2]
names_trs <- counts_transcripts[,c(1, 2, 15)]
# Selecting only counts
counts_genes <- counts_genes[, 3:14]
counts_transcripts <- counts_transcripts[, 3:14]

# Creating a DESeq2 dataset for GENES
ddsg <- DESeqDataSetFromMatrix(countData = round(counts_genes),
                              colData = data.frame(
                                condition = rep(c(
                                  "Control_AB1", "Control_AE21", "METTL1_KO_AB16", "METTL1_KO_AE4"
                                  ), each = 3)),
                              design = ~condition)
# Preliminary filtering, leaivng only entries that passed a minimum threshold
ddsg <- ddsg[ rowSums(counts(ddsg)) > 15, ]
ddsg <- ddsg[ rowMax(counts(ddsg)) > 10, ]
dim(counts(ddsg)) # 18701
sum(is.na(counts(ddsg)))

# Creating a DESeq2 dataset for TRANSCRIPTS
ddst <- DESeqDataSetFromMatrix(countData = round(counts_transcripts),
                               colData = data.frame(condition = rep(c(
                                 "Control_AB1", "Control_AE21", "METTL1_KO_AB16", "METTL1_KO_AE4"),
                                  each = 3)),
                               design = ~condition)
# Filtering
ddst <- ddst[ rowSums(counts(ddst)) > 15, ]
ddst <- ddst[ rowMax(counts(ddst)) > 10, ]
dim(counts(ddst))  # 88075
sum(is.na(counts(ddst)))

######################## DESeq2 model fitting, GENES ###########################

ddsg <- DESeq(ddsg)
head(results(ddsg))
summary(results(ddsg))
dim(assay(ddsg))

# Selecting KO-Ctrl comparisons
backgr_g <- results(ddsg, contrast = c("condition", "Control_AB1", "Control_AE21"))
res1_g <- results(ddsg, contrast = c("condition", "METTL1_KO_AE4", "Control_AB1"))
res2_g <- results(ddsg, contrast = c("condition", "METTL1_KO_AB16", "Control_AB1"))
res3_g <- results(ddsg, contrast = c("condition", "METTL1_KO_AE4", "Control_AE21"))
res4_g <- results(ddsg, contrast = c("condition", "METTL1_KO_AB16", "Control_AE21"))

##################### DESeq2 model fitting, TRANSCRIPTS ########################

ddst <- DESeq(ddst)
head(results(ddst))
summary(results(ddst))
dim(assay(ddst))

# Selecting KO-Ctrl comparisons
backgr_t <- results(ddst, contrast = c("condition", "Control_AB1", "Control_AE21"))
res1_t <- results(ddst, contrast = c("condition", "METTL1_KO_AE4", "Control_AB1"))
res2_t <- results(ddst, contrast = c("condition", "METTL1_KO_AB16", "Control_AB1"))
res3_t <- results(ddst, contrast = c("condition", "METTL1_KO_AE4", "Control_AE21"))
res4_t <- results(ddst, contrast = c("condition", "METTL1_KO_AB16", "Control_AE21"))

################# NA correction, adding back gene names ########################

# Function for NA in p-val handling, based on https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-are-some-p-values-set-to-na
pval_na_correction <- function(data){
  data <- as.data.frame(data)
  data$pvalue <- ifelse(is.na(data$pvalue), 1, data$pvalue)
  data$padj <- ifelse(is.na(data$padj), 1, data$padj)
  return(data)
}

# Adding the gene names
add_gene_names <- function(data){
  data$ids <- rownames(data)
  if(substr(rownames(data)[1], 4, 4) == "G") {
    data <- merge(data, names_genes, by.x = "ids", by.y = "gene_id")
  } else {
    data <- merge(data, names_trs, by.x = "ids", by.y = "tx")
  }
  rownames(data) <- data$ids
  data$ids <- NULL
  return(data)
}

data_frame_names <- c("backgr_t", "res1_t", "res2_t", "res3_t", "res4_t",
                      "backgr_g", "res1_g", "res2_g", "res3_g", "res4_g")
data_frame_list <- mget(data_frame_names)
corrected_data_frames <- lapply(data_frame_list, pval_na_correction)
corrected_data_frames2 <- lapply(corrected_data_frames, add_gene_names)

for (i in seq_along(data_frame_list)) {
  assign(data_frame_names[i], corrected_data_frames2[[i]])
}

# Checking the results 
sum(is.na(backgr_t))
head(backgr_t)
dim(backgr_g)

################# Visualization: MA plots, heatmaps for samples ################

# GENES:
sampleDists <- dist(t(assay(ddsg)))
sampleDistMatrix <- as.matrix(sampleDists)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         main = "Heatmap of euclidean distances between samples\n GENES")

# TRANSCRIPTS:
sampleDists <- dist(t(assay(ddst)))
sampleDistMatrix <- as.matrix( sampleDists )
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         main = "Heatmap of euclidean distances between samples\n TRANSCRIPTS")

# MA plots
plotMA(results(ddst), main = "MA Plot for transcripts", ylim=c(-4,4))
plotMA(results(ddsg), main = "MA Plot for genes", ylim=c(-4,4))

######################## Saving raw results, logFC sorted ######################

output_dir <- "./deseq2_results/"

# Function to export results to a csv file, sorted by absolute logFC
export_all_results_sorted <- function(results, tr_or_gene, comparison_name) {
  if(tr_or_gene == "gene") {
    type_dir <- 'genes'
  } else {type_dir <- 'transcripts'}
  all_genes_sorted <- results[order(abs(results$log2FoldChange), decreasing = TRUE), ]
  output_file <- file.path(output_dir, type_dir, paste0(comparison_name, "_deseq2_logFC_sorted.csv"))
  write.csv(all_genes_sorted, file = output_file, row.names = TRUE)
  cat("All results for", comparison_name, "exported to", output_file, "\n")
}

celllines_pairs <- list('OE33_METTL1_KO_AE4_vs_OE33_Control_AB1', 
                        'OE33_METTL1_KO_AB16_vs_OE33_Control_AB1', 
                        'OE33_METTL1_KO_AE4_vs_OE33_Control_AE21',
                        'OE33_METTL1_KO_AB16_vs_OE33_Control_AE21',
                        'OE33_Control_AB1_vs_OE33_Control_AE21')

deseq_raw_res_list_g <- list(res1_g, res2_g, res3_g, res4_g, backgr_g)
deseq_raw_res_list_t <- list(res1_t, res2_t, res3_t, res4_t, backgr_t)

# GENES:
for (i in seq_along(deseq_raw_res_list_g)) {
  export_all_results_sorted(deseq_raw_res_list_g[[i]], 'gene', celllines_pairs[[i]])
}

# TRANSCRIPTS:
for (i in seq_along(deseq_raw_res_list_t)) {
  export_all_results_sorted(deseq_raw_res_list_t[[i]], 'tr', celllines_pairs[[i]])
}

#################################### Filtering #################################

# Function to filter results by p-value and logFC and sort by absolute logFC
filter_and_sort <- function(data_init, p_val_threshold = 0.05, logFC_threshold = 0.5) {
  # for the DESeq2 analysis, we delete mitochondrial genes on this stage!
  data <- subset(data_init, !grepl("^MT-", gene_name))
  data <- data[data$pvalue < p_val_threshold, ]
  data <- data[abs(data$log2FoldChange) > logFC_threshold, ]
  data <- data[order(abs(data$log2FoldChange), decreasing = TRUE), ]
  print(paste('Size of filtered data: ', dim(data)))
  return(data)
}

deseq_filtered_g <- lapply(deseq_raw_res_list_g, filter_and_sort)
deseq_filtered_t <- lapply(deseq_raw_res_list_t, filter_and_sort)

# Saving filtered results for each comparison
type_dir <- 'transcripts'
for (i in seq_along(deseq_filtered_t)) {
  write.csv(deseq_filtered_t[[i]],
            file = file.path(output_dir, type_dir, paste0("Trs", i, "_filt_", celllines_pairs[[i]], "_deseq2.csv")),
            row.names = TRUE)
}
type_dir <- 'genes'
for (i in seq_along(deseq_filtered_g)) {
  write.csv(deseq_filtered_g[[i]],
            file = file.path(output_dir, type_dir, paste0("Genes", i, "_filt_", celllines_pairs[[i]], "_deseq2.csv")),
            row.names = TRUE)
}

################################# Saving names #################################

filter_and_names <- function(data_init, gene_names_for_tr = FALSE) {
  if(!gene_names_for_tr){
    if(substr(rownames(data_init)[1], 4, 4) == "G") {
      return(data_init$gene_name)
    } else {
      return(rownames(data_init))
    }
  } else {
    return(unique(data_init$gene_name))
  }
}

deseq_filtered_g_names <- lapply(deseq_filtered_g, filter_and_names)
deseq_filtered_t_names <- lapply(deseq_filtered_t, filter_and_names)

# GENES:
for (i in seq_along(deseq_filtered_g_names)) {
  write.table(deseq_filtered_g_names[[i]],
              file = paste0(output_dir,"genes/names/genes_",celllines_pairs[[i]],"_logfc05_pval005.txt"),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}
# TRANSCRIPTS:
for (i in seq_along(deseq_filtered_t_names)) {
  write.table(deseq_filtered_t_names[[i]],
              file = paste0(output_dir,"transcripts/names/trs_",celllines_pairs[[i]],"_logfc05_pval005.txt"),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}
