######################### Packages loading, setting up #########################

library(limma)
library(edgeR)
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

########################## Reading the data ####################################

# GENES
counts_genes <- read.delim("./2023_NovoGene_RNAseq_Report/OE33_salmon.merged.gene_counts.tsv")
# TRANSCRIPTS:
counts_transcripts <- read.delim("./2023_NovoGene_RNAseq_Report/OE33_salmon.merged.transcript_counts.tsv")

dim(counts_genes) # 62703
dim(counts_transcripts) # 252301
head(counts_genes)
head(counts_transcripts)
# Unique gene names - some gene names are still in ENSG form!
length(unique(counts_genes$gene_id))
# dim(subset(counts_genes, !grepl("^ENSG", gene_name)))

# Adding gene name column to the transcripts data
counts_transcripts <- merge(counts_transcripts, counts_genes[,1:2],
                            by.x = "gene_id", by.y = "gene_id")

# GENES: rownames are gene ids
rownames(counts_genes) <- counts_genes$gene_id
# TRANSCRIPTS: rownames are trs ids
rownames(counts_transcripts) <- counts_transcripts$tx

# Here I used the data with and without merged control reference from another 
# study available at https://www.ebi.ac.uk/ena/browser/view/SRX8604051
# Depending on what data was used, different indexes were selected.

# FOR DATA W/O REFERENCE
counts_genesL <- DGEList(counts=counts_genes[,3:14], genes=counts_genes[,1:2])
counts_transcriptsL <- DGEList(counts=counts_transcripts[,3:14], genes=counts_transcripts[,c(1:2, 15)])
# FOR DATA WITH REFERENCE
counts_geneL <- DGEList(counts=counts_gene[,c(2, 3, 4, seq(6,17))], genes=counts_gene[,c(1, 5)])
head(counts_genesL)
head(counts_transcriptsL)
counts_geneL$samples

################################# MT genes #####################################

# On the first step of data pre-processing we will remove mitochondrial genes
# since these genes require separated analysis. 

# Saving mt genes to the separated datasets
mt_transcripts <- subset(counts_transcriptsL$genes, grepl("^MT-", gene_name))
dim(mt_transcripts) # 37 trs
mt_genes <- subset(counts_genesL$genes, grepl("^MT-", gene_name))
dim(mt_genes) # 37 genes 

# Filtering out these genes from the main datasets
counts_transcriptsL <- NULL
counts_genesL <- NULL

genes_to_keep <- !rownames(counts_genesL$counts) %in% mt_genes$gene_id
counts_genesL$counts <- counts_genesL$counts[genes_to_keep, ]
counts_genesL$genes <- counts_genesL$genes[genes_to_keep, ]

trs_to_keep <- !rownames(counts_transcriptsL$counts) %in% mt_transcripts$tx
counts_transcriptsL$counts <- counts_transcriptsL$counts[genes_to_keep, ]
counts_transcriptsL$genes <- counts_transcriptsL$genes[genes_to_keep, ]

############################# Further filtering ################################

# Ordering genes
o <- order(rowSums(counts_genesL$counts), decreasing=TRUE)
counts_genesL <- counts_genesL[o,]
# Ordering transcripts
o <- order(rowSums(counts_transcriptsL$counts), decreasing=TRUE)
counts_transcriptsL <- counts_transcriptsL[o,]

# Removing duplicates
# GENES:
d <- duplicated(counts_genesL$genes$gene_name)
counts_genesL <- counts_genesL[!d,]
nrow(counts_genesL)  # 61248
# TRANSCRIPTS:
dt <- duplicated(counts_transcriptsL$genes$tx)
counts_transcriptsL <- counts_transcriptsL[!d,]
nrow(counts_transcriptsL)  # 246481

####################### Plotting raw expression values #########################

plot_expression <- function(de_data, tit){
  # Data preparation
  melted_counts <- melt(de_data$counts)
  merged_data <- merge(melted_counts, de_data$genes, by.x = "Var1",
                       by.y = "row.names", all.x = TRUE)
  #head(merged_data)
  colnames(merged_data) <- c("Gene_ID", "Sample", "Count", "Gene_id_old", "Gene_Name")
  
  # Calculating mean values for each sample
  mean_data <- merged_data %>%
    group_by(Sample) %>%
    summarize(mean_count = mean(Count, na.rm = TRUE))
  
  # Plotting a expression distribution plot for each sample
  expression_plots <- ggplot(merged_data, aes(x = log(Count))) +
    geom_histogram(bins = 40) +
    geom_vline(xintercept = log(10), linetype = "dashed",
               color = "red", size = 0.5, alpha = 0.5) +
    geom_text(data = mean_data, aes(x = Inf, y = Inf, label = sprintf("Mean: %.2f", mean_count)),
              hjust = 1, vjust = 1, size = 3, color = "black") +
    facet_wrap(~Sample, scales = "free_y", ncol=2) +
    labs(title = tit,
         x = "Log(expression)",
         y = "Frequency") +
    theme_minimal()
  
  print(expression_plots)
}

plot_expression(counts_geneL, "Log expression distribution before filtering")

############################ Creating the design matrix ########################

# For data with reference sample
names <- sub("_R.*", "", colnames(counts_gene)[6:17])
names <- c('ref', 'ref', 'ref', names)
counts_geneL$samples$group <- names
# For data w/o reference
names <- sub("_R.*", "", colnames(counts_genes)[3:14])
counts_genesL$samples$group <- names
counts_transcriptsL$samples$group <- names

design_matrix <- model.matrix(~0 + names)

######################### Filtering low expression genes #######################

# GENES:
keep <- filterByExpr(counts_genesL, design = design_matrix)
counts_genesL <- counts_genesL[keep, , keep.lib.sizes=FALSE]
nrow(counts_genesL) # 17322
# TRANSCRIPTS:
keep <- filterByExpr(counts_transcriptsL, design = design_matrix)
counts_transcriptsL <- counts_transcriptsL[keep, , keep.lib.sizes=FALSE]
nrow(counts_transcriptsL) # 68722

plot_expression(counts_geneL, "Log expression distribution after filtering")

# Recomputing library size
counts_genesL$samples$lib.size <- colSums(counts_genesL$counts)
counts_transcriptsL$samples$lib.size <- colSums(counts_transcriptsL$counts)

################################# Normalization ################################

# TMM normalization to account for compositional differences between the libraries
counts_genesL <- calcNormFactors(counts_genesL)
counts_transcriptsL <- calcNormFactors(counts_transcriptsL)
counts_genesL$samples

short_names <- c(#'ref1', 'ref2', 'ref3', <- with ref
                 "cAB1", "cAB2", "cAB3",
                 "cAE1", "cAE2", "cAE3",
                 "mAB1", "mAB2", "mAB3",
                 "mAE1", "mAE2", "mAE3")

plotMDS(counts_genesL, labels = short_names)
  
############################### Defining contrasts #############################

contrast_matrix <- makeContrasts(
  OE33_METTL1_KO_AE4_vs_OE33_Control_AB1 = namesOE33_METTL1_KO_AE4 - namesOE33_Control_AB1,
  OE33_METTL1_KO_AB16_vs_OE33_Control_AB1 = namesOE33_METTL1_KO_AB16 - namesOE33_Control_AB1,
  OE33_METTL1_KO_AE4_vs_OE33_Control_AE21 = namesOE33_METTL1_KO_AE4 - namesOE33_Control_AE21,
  OE33_METTL1_KO_AB16_vs_OE33_Control_AE21 = namesOE33_METTL1_KO_AB16 - namesOE33_Control_AE21,
  # the contrast between two controls:
  OE33_Control_AB1_vs_OE33_Control_AE21 = namesOE33_Control_AB1 - namesOE33_Control_AE21,
  #ref_vs_OE33_Control_AB1 = namesref - namesOE33_Control_AB1,   <- with ref
  #ref_vs_OE33_Control_AE21 = namesref - namesOE33_Control_AE21,   <- with ref
  levels = design_matrix
)

# Estimating dispersion
# GENES:
counts_genesL <- estimateGLMCommonDisp(counts_genesL)
counts_genesL <- estimateGLMTrendedDisp(counts_genesL)
counts_genesL <- estimateGLMTagwiseDisp(counts_genesL)
# TRANSCRIPTS:
counts_transcriptsL <- estimateGLMCommonDisp(counts_transcriptsL)
counts_transcriptsL <- estimateGLMTrendedDisp(counts_transcriptsL)
counts_transcriptsL <- estimateGLMTagwiseDisp(counts_transcriptsL)

# Fit the negative binomial generalized log-linear model
fit_genes <- glmFit(counts_genesL, design_matrix)
fit_tr <- glmFit(counts_transcriptsL, design_matrix)

############################### Contrast analysis ##############################

num_contrasts <- 5  # 7 for data with reference
contrast_list_genes <- vector("list", num_contrasts) 
contrast_list_trs <- vector("list", num_contrasts)

# Saving computed values for each comparison
for (i in 1:num_contrasts) {
  contrast_list_genes[[i]] <- glmLRT(fit_genes, contrast = contrast_matrix[, i])
  contrast_list_trs[[i]] <- glmLRT(fit_tr, contrast = contrast_matrix[, i])
}

# Adjusting p-values for each contrast
for (i in 1:num_contrasts) {
  contrast_list_genes[[i]]$table$FDR <- p.adjust(contrast_list_genes[[i]]$table$PValue, method = "BH")
  contrast_list_trs[[i]]$table$FDR <- p.adjust(contrast_list_trs[[i]]$table$PValue, method = "BH")
}

############################# Saving results to files ##########################

output_dir <- "./edgeR_results/"

# Function for saving results sorted by logFC to the csv file
# (not filtered by p-value and logFC!)
export_all_results_sorted <- function(results, comparison_name) {
  all_genes <- as.data.frame(results)
  all_genes_sorted <- all_genes[order(abs(all_genes$logFC), decreasing = TRUE), ]
  output_file <- file.path(output_dir, paste0(comparison_name, "_edgeR_logFC_sorted.csv"))
  write.csv(all_genes_sorted, file = output_file, row.names = TRUE)
  cat("All results for", comparison_name, "exported to", output_file, "\n")
}

celllines_pairs <- list('OE33_METTL1_KO_AE4_vs_OE33_Control_AB1', 
                        'OE33_METTL1_KO_AB16_vs_OE33_Control_AB1', 
                        'OE33_METTL1_KO_AE4_vs_OE33_Control_AE21',
                        'OE33_METTL1_KO_AB16_vs_OE33_Control_AE21',
                        'OE33_Control_AB1_vs_OE33_Control_AE21')

# GENES:
for (i in seq_along(contrast_list_genes)) {
  export_all_results_sorted(contrast_list_genes[[i]],
                            paste0(output_dir, "Genes_", celllines_pairs[[i]]))
}

# TRANSCRIPTS:
for (i in seq_along(contrast_list_trs)) {
  export_all_results_sorted(contrast_list_trs[[i]],
                            paste0(output_dir, "Trs_", celllines_pairs[[i]]))
}

################ Converting results to dataframes, filtering ###################

# Function to filter results by p-value and logFC and sort by absolute logFC
filter_and_sort <- function(data_init, p_val_threshold = 0.05, logFC_threshold = 0.5) {
  data <- as.data.frame(data_init$table)
  data$gene_name <- data_init$genes$gene_name
  data <- data[data$PValue < p_val_threshold, ]
  data <- data[abs(data$logFC) > logFC_threshold, ]
  data <- data[order(abs(data$logFC), decreasing = TRUE), ]
  print(paste('Size of filtered data: ', dim(data)))
  return(data)
}

filtered_contrasts_g <- lapply(contrast_list_genes, filter_and_sort)
filtered_contrasts_t <- lapply(contrast_list_trs, filter_and_sort)

# Saving filtered results for each comparison
for (i in seq_along(filtered_contrasts_t)) {
  write.csv(filtered_contrasts_t[[i]],
            file = paste0(output_dir, "Trs", i, "_filt_", celllines_pairs[[i]], "_edgeR.csv"),
            row.names = TRUE)
}

for (i in seq_along(filtered_contrasts_g)) {
  write.csv(filtered_contrasts_g[[i]],
            file = paste0(output_dir, "Genes", i, "_filt_", celllines_pairs[[i]], "_edgeR.csv"),
            row.names = TRUE)
}

############################## Visualization: PCA ##############################

# Calculating scaling factors to convert raw library sizes into effective library sizes
norm_factors <- calcNormFactors(counts_geneL)
counts_normalized <- cpm(counts_geneL)

pca_result <- prcomp(t(counts_normalized))
pca_data <- as.data.frame(pca_result$x)
pca_data$group <- counts_geneL$samples$group

ggplot(pca_data, aes(x = PC1, y = PC2, color = group)) +
  geom_point() +
  labs(title = "PCA plot for METTL1 KO and Control OE33 samples")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

######################## Visualization: tSNE & UMAP ############################

### tSNE & UMAP on the raw data

# t-SNE on raw data
tsne_result_raw <- Rtsne(t(as.matrix(counts_normalized)), perplexity = 3)
# UMAP on raw data
umap_result_raw <- umap(t(as.matrix(counts_normalized)), n_neighbors = 5)

tsne_umap_data_raw <- data.frame(
  tSNE1 = tsne_result_raw$Y[, 1],
  tSNE1 = tsne_result_raw$Y[, 2],
  UMAP1 = umap_result_raw$layout[, 1],
  UMAP2 = umap_result_raw$layout[, 2],
  group = counts_geneL$samples$group
)

# Ploting
ggplot(tsne_umap_data_raw, aes(x = tSNE1, y = tSNE2, color = group)) +
  geom_point() +
  labs(title = "t-SNE for METTL1 KO and Control OE33 samples, raw") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(tsne_umap_data_raw, aes(x = UMAP1, y = UMAP2, color = group)) +
  geom_point() +
  labs(title = "UMAP for METTL1 KO and Control OE33 samples, raw") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

### tSNE & UMAP on the PCA-transformed data

# t-SNE
tsne_result <- Rtsne(as.matrix(pca_data[, -ncol(pca_data)]), perplexity = 3)
# UMAP
umap_result <- umap(as.matrix(pca_data[, -ncol(pca_data)]), n_neighbors = 3)

# Add t-SNE and UMAP results to PCA data
tsne_umap_data <- data.frame(
  TSNE1 = tsne_result$Y[, 1],
  TSNE2 = tsne_result$Y[, 2],
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  group = counts_geneL$samples$group
)

# Plotting
ggplot(tsne_umap_data, aes(x = TSNE1, y = TSNE2, color = group)) +
  geom_point() +
  labs(title = "t-SNE for METTL1 KO and Control OE33 samples, on PCA") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(tsne_umap_data, aes(x = UMAP1, y = UMAP2, color = group)) +
  geom_point() +
  labs(title = "UMAP for METTL1 KO and Control OE33 samples, on PCA") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

######################## Visualization: Volcano Plots ##########################

# Function to create and save volcano plot
create_and_save_volcano_plot <- function(contrast_result, filename) {
  logFC <- contrast_result$table$logFC
  logPValue <- -log10(contrast_result$table$PValue)
  gene_ids <- contrast_result$genes$gene_name
  colors <- ifelse(contrast_result$table$PValue < 0.05 & logFC > 2, "red",
                   ifelse(contrast_result$table$PValue < 0.05 & logFC < -2, "green", "black"))
  
  volcano_data <- data.frame(logFC, logPValue, gene_ids, colors)
  
  volcano_plot <- ggplot(volcano_data, aes(x = logFC, y = logPValue, color = colors)) +
    geom_point(alpha = 0.5) +
    scale_color_identity() +
    labs(title = gsub(".png", "", filename), x = "log2 Fold Change", y = "-log10(PValue)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +  # Center the title
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
    geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "gray") +
    geom_text_repel(data = subset(volcano_data, colors != "black"), aes(label = gene_ids), 
                    box.padding = 0.5, point.padding = 0.1, size = 3, color = "black")  # Repel gene names
  
  ggsave(filename, plot = volcano_plot, device = "png", width = 8, height = 6, units = "in")
}

output_dir_img <-'./images/volcano_plots/edgeR/'

# GENES:
for (i in seq_along(contrast_list_genes)) {
  create_and_save_volcano_plot(contrast_list_genes[[i]],
                            paste0(output_dir_img, "Genes_", celllines_pairs[[i]]))
}

# TRANSCRIPTS:
for (i in seq_along(contrast_list_trs)) {
  create_and_save_volcano_plot(contrast_list_trs[[i]],
                            paste0(output_dir_img, "Trs_", celllines_pairs[[i]]))
}

############################# Saving names of top genes ########################

# filtering more strictly
#filter_and_names <- function(data_init, logFC_trh = 1, FRD_trh = 0.2) {
  # Filtering by p-value and logFC
  #data <- data_init[(data_init$FDR < FRD_trh) & (abs(data_init$logFC) > logFC_trh), ]
  #print(dim(data))
  #return(data$gene_name)
#}

select_names <- function(data_init, gene_names_for_tr = FALSE) {
  if(!gene_names_for_tr){
    if(substr(rownames(data_init)[1], 4, 4) == "G") {
      return(data_init$gene_name) # returning gene names
    } else {
      return(rownames(data_init)) # returning transcript names
    }
  } else {
    return(unique(data_init$gene_name)) # returning gene names for trs as well
  }
}

filtered_contrasts_g_names <- lapply(filtered_contrasts_g, select_names)
filtered_contrasts_t_names <- lapply(filtered_contrasts_t, select_names)

# GENES:
for (i in seq_along(filtered_contrasts_g_names)) {
  write.table(filtered_contrasts_g_names[[i]],
              file = paste0(output_dir,"genes/names/genes_",celllines_pairs[[i]],"_logfc05_pval005.txt"),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}
# TRANSCRIPTS:
for (i in seq_along(filtered_contrasts_t_names)) {
  write.table(filtered_contrasts_t_names[[i]],
              file = paste0(output_dir,"trs/names/trs_",celllines_pairs[[i]],"_logfc05_pval005.txt"),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}
