###################### Packages loading, setting up ############################

library(dplyr)
library(stringr)
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
library(eulerr)

setwd("~/Downloads/OE33_project")

#################### Full RNA-seq & full proteomics intersection ###############

# Raw data intersections - selecting common names
full_rnaseq_prots_names_t <- intersect(contrast_list_trs[[1]]$genes$gene_name,
                                       intens_for_kp$symbol) # 3329
full_rnaseq_prots_names_g <- intersect(contrast_list_genes[[1]]$genes$gene_name,
                                       intens_for_kp$symbol) # 3331
length(full_rnaseq_prots_names_g)

# Intersections with the comparisons filtered by p-val<0.05 and lfc>0.5, 
# using data-intersection of edgeR and DESeq2

gene_int <- c('ae4_ab1_intg', 'ab16_ab1_intg', 'ae4_ae21_intg',
                        'ab16_ae21_intg', 'ab1_ae21_intg')
trs_int <- c('ae4_ab1_intt', 'ab16_ab1_intt', 'ae4_ae21_intt',
                              'ab16_ae21_intt', 'ab1_ae21_intt')

rna_prot_inter_genes <- list()
rna_prot_inter_trs <- list()

for (i in 1:5) {
  rna_prot_inter_genes[[gene_int[i]]] <- intersect(ints_genes_full_list[[gene_int[i]]]$gene_name,
                                                   intens_for_kp$symbol)
  rna_prot_inter_trs[[trs_int[i]]] <- intersect(ints_trs_full_list[[trs_int[i]]]$gene_name,
                                                intens_for_kp$symbol)
}

# Number of proteins in this intersection:
# with GENES:
# AE4-AB1: 886
# AB16-AB1: 571
# AB16-AE21: 760
# AE4-AE21: 1179
# with TRANSCRIPS:
# AE4-AB1: 427
# AB16-AB1: 279
# AB16-AE21: 349
# AE4-AE21: 504

# Intersection of the proteomics with the 'global intersection'
# (consistently up-/down-regulated genes/trs)

gi_g_pinters <-intersect(global_int_g_counts$gene_name, intens_for_kp$symbol) # 114
gi_t_pinters <-intersect(global_int_t_counts$gene_name, intens_for_kp$symbol) # 34

# In order to take the largest set of proteins for the network analysis,
# but at the same time consider only consistent proteins, we will take a union
# of the AE4-AB1, AB16-AB1, AB16-AE21, and AE4-AE21 lists. By this we are selecting 
# genes/trs that are 1) in proteomics 2) show p-val < 0.05 and |logFC| > 0.5 in at 
# least of one KO-Ctrl comparison.

genes_union <- union(union(ae4_ab1_gpinters, ab16_ab1_gpinters), union(ab16_ae21_gpinters, ae4_ae21_gpinters))
length(genes_union) # 1821
trs_union <- union(union(ae4_ab1_tpinters, ab16_ab1_tpinters), union(ab16_ae21_tpinters, ae4_ae21_tpinters))
length(trs_union) # 946

# Some genes from these unions, however, could have different directions of logFC
# in different comparisons. Here we will filter out these genes.

down_genes <- select_one_direction(genes_union)[['down']] # 598
up_genes <- select_one_direction(genes_union)[['up']] # 1167
down_trs <- select_one_direction(trs_union)[['down']] # 283
up_trs <- select_one_direction(trs_union)[['up']] # 633
genes_union <- append(down_genes, up_genes)
trs_union <- append(down_trs, up_trs)

# We will perform the same selection of entries with the same direction but now 
# for proteomics:

dfg <- data.frame(genes_union)
dft <- data.frame(trs_union)

# Adding the column indicating for each protein if the DE change for it in proteomics 
# is in the same direction for AB16-NC2 and AE4-NC2 comparisons.

dfg$same_direction <- sapply(dfg$genes_union, function(name){
  sign(subset(res_AE4, res_AE4$protein==name)$logFC) == sign(subset(res_AB16, res_AB16$protein==name)$logFC)
})
dfg$same_direction_imp <- sapply(dfg$genes_union, function(name){
  sign(subset(res_AE4_imp, res_AE4_imp$protein==name)$logFC) == sign(subset(res_AB16_imp, res_AB16_imp$protein==name)$logFC)
})
dft$same_direction <- sapply(dft$trs_union, function(name){
  sign(subset(res_AE4, res_AE4$protein==name)$logFC) == sign(subset(res_AB16, res_AB16$protein==name)$logFC)
})
dft$same_direction_imp <- sapply(dft$trs_union, function(name){
  sign(subset(res_AE4_imp, res_AE4_imp$protein==name)$logFC) == sign(subset(res_AB16_imp, res_AB16_imp$protein==name)$logFC)
})

# Saving only proteins with the common direction
df_fin_t <- filter(dft, same_direction_imp==TRUE)
dim(df_fin_t) # 557
df_fin_g <- filter(dfg, same_direction_imp==TRUE)
dim(df_fin_g) # 1055

write.table(df_fin_t$trs_union,
            file = "early_prot_inters_trs_names.txt",
            row.names = FALSE, quote = F, col.names = F)
write.table(df_fin_g$genes_union,
            file = "early_prot_inters_genes_names.txt",
            row.names = FALSE, quote = F, col.names = F)

# Adding mean logFC to this data for the SPIA analysis

source1 <- res_AB16_imp
source2 <- res_AE4_imp
df_fin_t$mean_lfc <- sapply(df_fin_t$trs_union, function(name) {
  mean(c(subset(source1, source1$protein==name)$logFC,
         subset(source2, source2$protein==name)$logFC))
})
df_fin_g$mean_lfc <- sapply(df_fin_g$genes_union, function(name) {
  mean(c(subset(source1, source1$protein==name)$logFC,
         subset(source2, source2$protein==name)$logFC))
})

# Intersection between transcripts, proteomics, and genes 

# There are 479 proteins appearing in both proteomics-trs and proteomics-genes lists
three_inters <- data.frame()
names_in_trs_and_genes <- intersect(df_fin_t$trs_union, df_fin_g$genes_union)

for(i in 1:nrow(df_fin_t)) {
  # we add the protein if it is 1) in both proteomics-trs and proteomics-genes lists
  # and 2) have the same DE change direction there
  if((df_fin_t$trs_union[i] %in% names_in_trs_and_genes) && (sign(df_fin_t$mean_lfc[i]) == sign(df_fin_g$mean_lfc[i]))){
    row_to_add <- data.frame(gene_name = df_fin_t$trs_union[i], mean_lfc_tr = df_fin_t$mean_lfc[i], mean_lfc_g = df_fin_g$mean_lfc[i])
    three_inters <- rbind(three_inters, row_to_add)
  }
}

# There are 240 genes in the proteomics & rna_seq_genes & rna_seq_transcripts
# intersection that are changing in the same direction when comparing between 
# METTL1-KO and control. 

write.table(three_inters$gene_name,
            file = "early_prot_3inters_names.txt",
            row.names = FALSE, quote = F, col.names = F)


##### Which key genes are also in proteomics? heatmap (w/o parental samples) #####

key_prots_proteomics <- intersect(key_prots, unique(intens$symbol)) # 28 proteins
key_intens <- subset(intens, symbol %in% key_prots)
key_intens$accession <- NULL
dim(key_intens)
head(key_intens)

# we want to throw out the parental wt
key_intens <- key_intens[,c(1, 11:37)]
key_intens <- subset(key_intens, symbol != "TYMS")

names(key_intens) <- c("gene_name",
                       # "Parental_WT_1.1","Parental_WT_1.2","Parental_WT_1.3",
                       # "Parental_WT_2.1","Parental_WT_2.2","Parental_WT_2.3",
                       # "Parental_WT_3.1","Parental_WT_3.2","Parental_WT_3.3",
                       "AB16_METTL1_KO_1.1","AB16_METTL1_KO_1.2","AB16_METTL1_KO_1.3",
                       "AB16_METTL1_KO_2.1","AB16_METTL1_KO_2.2","AB16_METTL1_KO_2.3",
                       "AB16_METTL1_KO_3.1","AB16_METTL1_KO_3.2","AB16_METTL1_KO_3.3",
                       "AE4_METTL1_KO_1.1","AE4_METTL1_KO_1.2","AE4_METTL1_KO_1.3",
                       "AE4_METTL1_KO_2.1","AE4_METTL1_KO_2.2","AE4_METTL1_KO_2.3",
                       "AE4_METTL1_KO_3.1","AE4_METTL1_KO_3.2","AE4_METTL1_KO_3.3",
                       "NC2_WT_1.1","NC2_WT_1.2","NC2_WT_1.3",
                       "NC2_WT_2.1","NC2_WT_2.2","NC2_WT_2.3",
                       "NC2_WT_3.1","NC2_WT_3.2","NC2_WT_3.3")

# for col clusters
samples <- data.frame(sample = as.character(names(key_intens[,2:28])), Sample_type="METTL1_KO") %>%
  column_to_rownames("sample")
samples$Sample_type <- c(rep("AB16_METTL1_KO", 9),
                         rep("AE4_METTL1_KO", 9), rep("NC2_WT", 9))

# for row clusters
up_down_genes <- data.frame(row = as.character(rownames(key_intens)),
                            gene_name = as.character(key_intens$gene_name))
rownames(up_down_genes) <- up_down_genes$gene_name
up_down_genes$Regulation <- ifelse(up_down_genes$gene_name %in% down_genes, "Down", "Up")

o <- order(up_down_genes$Regulation)
up_down_genes <- up_down_genes[o, ]
rownames(key_intens) <- key_intens$gene_name
key_intens <- key_intens[, 2:28]
up_down_genes$row <- NULL
up_down_genes$gene_name <- NULL

hm_gg <- as.ggplot(pheatmap(as.matrix(log(key_intens+1)), scale="row",
                            cluster_rows = T,
                            cluster_cols = F,
                            fontsize_col = 8,
                            fontsize_row = 8,
                            angle_col = 90,
                            cellheight = 10,
                            treeheight_col = 0,
                            treeheight_row = 4,
                            annotation_row = up_down_genes,
                            annotation_col= samples,
                            annotation_colors =list(
                              Sample_type=c(AB16_METTL1_KO="darkgoldenrod1",
                                            AE4_METTL1_KO="coral3", NC2_WT="bisque4"),
                              Regulation=c(Up="chocolate1", Down="cornflowerblue")),
                            color=colorRampPalette(c("navy","white", "red"))(50), 
                            show_rownames = T, annotation_names_col = F, annotation_names_row = F))


######## Which key genes are also in proteomics? heatmap (only 3rd reps) ########

key_prots_proteomics <- intersect(key_prots, unique(intens$symbol)) # 28 proteins
key_intens <- subset(intens, symbol %in% key_prots)
key_intens$accession <- NULL
dim(key_intens)
head(key_intens)

# we want to LEAVE ONLY 3RD REPS AND throw out parental
key_intens <- key_intens[,c(1, 17:19, 26:28, 35:37)]
key_intens <- subset(key_intens, symbol != "TYMS")
key_intens <- subset(key_intens, symbol != "PKM")
# leaving only the consistently upregulated prots
key_intens <- key_intens[apply(key_intens[, 2:7], 1, min) > apply(key_intens[, 8:10], 1, max), ]

names(key_intens) <- c("gene_name","AB16_METTL1_KO_3.1","AB16_METTL1_KO_3.2",
                       "AB16_METTL1_KO_3.3","AE4_METTL1_KO_3.1","AE4_METTL1_KO_3.2",
                       "AE4_METTL1_KO_3.3","NC2_WT_3.1","NC2_WT_3.2","NC2_WT_3.3")

# for col clusters
samples <- data.frame(sample = as.character(names(key_intens[, 2:10])), Sample_type="METTL1_KO") %>%
  column_to_rownames("sample")
samples$Sample_type <- c(rep("AB16_METTL1_KO", 3),
                         rep("AE4_METTL1_KO", 3), rep("NC2_WT", 3))

# for row clusters
up_down_genes <- data.frame(row = as.character(rownames(key_intens)),
                            gene_name = as.character(key_intens$gene_name))
rownames(up_down_genes) <- up_down_genes$gene_name
up_down_genes$Regulation <- ifelse(up_down_genes$gene_name %in% down_genes, "Down", "Up")

o <- order(up_down_genes$Regulation)
up_down_genes <- up_down_genes[o, ]
rownames(key_intens) <- key_intens$gene_name
key_intens <- key_intens[, 2:10]
up_down_genes$row <- NULL
up_down_genes$gene_name <- NULL

hm_gg <- as.ggplot(pheatmap(as.matrix(log(key_intens+1)), scale="row",
                            cluster_rows = T,
                            cluster_cols = F,
                            fontsize_col = 8,
                            fontsize_row = 8,
                            angle_col = 90,
                            cellheight = 10,
                            treeheight_col = 0,
                            treeheight_row = 4,
                            annotation_row = up_down_genes,
                            annotation_col= samples,
                            annotation_colors =list(
                              Sample_type=c(AB16_METTL1_KO="darkgoldenrod1",
                                            AE4_METTL1_KO="coral3", NC2_WT="bisque4"),
                              Regulation=c(Up="chocolate1", Down="cornflowerblue")),
                            color=colorRampPalette(c("navy","white", "red"))(50), 
                            show_rownames = T, annotation_names_col = F, annotation_names_row = F))

