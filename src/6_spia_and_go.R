############################### SPIA, libraries ################################

library(limma)
library(edgeR)
library(dplyr)
library(reshape2)
library(tibble)
library(tidyverse)

library(ggplot2)
library(gplots)
library(ggrepel)
library(ggplotify)
library(pheatmap)
library(heatmaply)
library(RColorBrewer)
library(grDevices)
library(eulerr)

library(SPIA)
library(GO.db)
library(org.Hs.eg.db)

setwd("~/Downloads/OE33_project")

##################### SPIA: for all comparisons separately #####################

# Filtering genes by p-value and logFC
filt_de <- function(table_df, genes_df, output_file) {
  combined_df <- cbind(table_df, genes_df)
  filtered_genes <- combined_df[((combined_df$PValue < 0.05) & (abs(combined_df$logFC) > 2)), ]
  sorted_genes <- filtered_genes[order(filtered_genes$logFC, decreasing = TRUE), ]
  return(sorted_genes)
}

edgeR_c1_filt_de <- filt_de(contrast_list_genes[[1]]$table, contrast_list_genes[[1]]$genes)
edgeR_c2_filt_de <- filt_de(contrast_list_genes[[2]]$table, contrast_list_genes[[2]]$genes)
edgeR_c3_filt_de <- filt_de(contrast_list_genes[[3]]$table, contrast_list_genes[[3]]$genes)
edgeR_c4_filt_de <- filt_de(contrast_list_genes[[4]]$table, contrast_list_genes[[4]]$genes)

# Below I will save filtered datasets for the each KO-Ctrl comparison to the 
# separated variable to compare it to the whole list of genes. Everywhere the 
# ENTREZ gene ID is used.

# Comparison 1: AE4-AB1
edgeR_c1_filt_de$ENTREZ <- mapIds(org.Hs.eg.db, keys = edgeR_c1_filt_de$gene_name,
                                  column = "ENTREZID", keytype = "SYMBOL")
edgeR_c1_filt_de <- edgeR_c1_filt_de[!is.na(edgeR_c1_filt_de$ENTREZ),]
DE_c1_e <- edgeR_c1_filt_de$logFC
names(DE_c1_e) <- as.vector(edgeR_c1_filt_de$ENTREZ)

# Comparison 2: AB16-AB1
edgeR_c2_filt_de$ENTREZ <- mapIds(org.Hs.eg.db, keys = edgeR_c2_filt_de$gene_name,
                                  column = "ENTREZID", keytype = "SYMBOL")
edgeR_c2_filt_de <- edgeR_c1_filt_de[!is.na(edgeR_c2_filt_de$ENTREZ),]
DE_c2_e <- edgeR_c2_filt_de$logFC
names(DE_c2_e) <- as.vector(edgeR_c2_filt_de$ENTREZ)

# Comparison 3: AE4-AE21
edgeR_c3_filt_de$ENTREZ <- mapIds(org.Hs.eg.db, keys = edgeR_c3_filt_de$gene_name,
                                  column = "ENTREZID", keytype = "SYMBOL")
edgeR_c3_filt_de <- edgeR_c3_filt_de[!is.na(edgeR_c3_filt_de$ENTREZ),]
DE_c3_e <- edgeR_c3_filt_de$logFC
names(DE_c3_e) <- as.vector(edgeR_c3_filt_de$ENTREZ)

# Comparison 4: AB16-AE21
edgeR_c4_filt_de$ENTREZ <- mapIds(org.Hs.eg.db, keys = edgeR_c4_filt_de$gene_name,
                                  column = "ENTREZID", keytype = "SYMBOL")
edgeR_c4_filt_de <- edgeR_c1_filt_de[!is.na(edgeR_c4_filt_de$ENTREZ),]
DE_c4_e <- edgeR_c4_filt_de$logFC
names(DE_c4_e) <- as.vector(edgeR_c4_filt_de$ENTREZ)

# ALL genes, the background to compare to
counts_gene$ENTREZ <- mapIds(org.Hs.eg.db, keys = counts_gene$gene_name,
                             column = "ENTREZID", keytype = "SYMBOL")
counts_gene_na <- counts_gene[!is.na(counts_gene$ENTREZ),]
ALL <- counts_gene_na$ENTREZ

# Performing a pathway SPIA analysis based on combined evidence
# (from help: use nB=2000 or more for more accurate results)

# Comparison 1:
res <- spia(de=DE_c1_e, all=ALL, organism="hsa",
            nB=2000,plots=FALSE,beta=NULL,combine="fisher",verbose=FALSE)
res[1:10, c(1, 11)]
write.csv(res, file = "spia_c1_res.csv", row.names = TRUE)

# Comparison 2:
res2 <- spia(de=DE_c2_e, all=ALL, organism="hsa",
             nB=2000,plots=FALSE,beta=NULL,combine="fisher",verbose=FALSE)
res2[1:10, c(1, 11)]
write.csv(res2, file = "spia_c2_res.csv", row.names = TRUE)

# Comparison 3:
res3 <- spia(de=DE_c3_e, all=ALL, organism="hsa",
             nB=2000,plots=FALSE,beta=NULL,combine="fisher",verbose=FALSE)
res3[1:10, c(1, 11)]
write.csv(res3, file = "spia_c3_res.csv", row.names = TRUE)

# Comparison 4:
res4 <- spia(de=DE_c4_e, all=ALL, organism="hsa",
             nB=2000,plots=FALSE,beta=NULL,combine="fisher",verbose=FALSE)
res4[1:10, c(1, 11)]
write.csv(res4, file = "spia_c4_res.csv", row.names = TRUE)

# Visualize results (optional)
#plotPathwayGraph(pathway_data$graph, highlightNodes = spia_result$id)


####################### SPIA: for all DE genes/trs #############################

###### GENES:

# All genes to compare to (background) - just preliminary filtering
counts_df_g <- as.data.frame(fit_genes$counts) # 17322
counts_df_g$gene_id <- as.character(fit_genes$genes$gene_id)
dim(counts_df_g) # 17322

# Adding mean logFC column to the full dataset
counts_df_g$mean_lfc <- sapply(counts_df_g$gene_id, function(name){
  mean(c(subset(contrast1g$table, rownames(contrast1g$table)==name)$logFC,
         subset(contrast2g$table, rownames(contrast2g$table)==name)$logFC,
         subset(contrast3g$table, rownames(contrast3g$table)==name)$logFC,
         subset(contrast4g$table, rownames(contrast4g$table)==name)$logFC))
})
counts_df_g$ENTREZ <- mapIds(org.Hs.eg.db, keys = counts_df_g$gene_id,
                             column = "ENTREZID", keytype = "ENSEMBL")
ALL <- unique(counts_df_g[!is.na(counts_df_g$ENTREZ),]$ENTREZ)
dim(counts_df_g)  # 14873

# Selecting DE genes: for union(DE genes from all KO-Ctrl comps) - Ctrl-Ctrl
all_de_names <- setdiff(union(union(rownames(comp1g), rownames(comp2g)),
                              union(rownames(comp3g), rownames(comp4g))), rownames(comp5g)) # 7749
all_de_df <- as.data.frame(all_de_names)
head(all_de_df)
all_de_df$mean_lfc <- sapply(all_de_df$all_de_names, function(name) {
  ifelse(name %in% counts_df_g$gene_id, 
         counts_df_g$mean_lfc[which(counts_df_g$gene_id == name)], 
         0)
})
# Adding the entrez id
all_de_df$ENTREZ <- mapIds(org.Hs.eg.db, keys = all_de_df$all_de_names,
                           column = "ENTREZID", keytype = "ENSEMBL")

# Leaving excising entrez names + highest value among the duplicated value 
all_de_df <- all_de_df[!is.na(all_de_df$ENTREZ),]
all_de_df$mean_lfc <- as.numeric(all_de_df$mean_lfc)
all_de_df <- all_de_df %>%
  group_by(ENTREZ) %>%
  slice_max(order_by = mean_lfc, n = 1) %>%
  ungroup()
all_de_spia_list <- all_de_df$mean_lfc
names(all_de_spia_list) <- as.vector(all_de_df$ENTREZ)

# Performing SPIA DE-vs-all
res_de <- spia(de=all_de_spia_list, all=ALL, organism="hsa",
               nB=500,plots=FALSE,beta=NULL,combine="fisher",verbose=FALSE)
res_de[1:10, c(1, 11)]
write.csv(res_de, file = "../spia/spia_res_rnaseq_fullDE_genes.csv", row.names = TRUE)


###### TRANSCRIPTS:

# All trs to compare to (background) - just preliminary filtering
counts_df_t <- as.data.frame(fit_tr$counts) # 68722
counts_df_t$tx_id <- as.character(fit_tr$genes$tx)
dim(counts_df_t) # 68722

# Adding mean logFC column to the full dataset
counts_df_t$mean_lfc <- sapply(counts_df_t$tx_id, function(name){
  mean(c(subset(contrast1t$table, rownames(contrast1t$table)==name)$logFC,
         subset(contrast2t$table, rownames(contrast2t$table)==name)$logFC,
         subset(contrast3t$table, rownames(contrast3t$table)==name)$logFC,
         subset(contrast4t$table, rownames(contrast4t$table)==name)$logFC))
})
counts_df_t$ENTREZ <- mapIds(org.Hs.eg.db, keys = counts_df_t$tx_id,
                             column = "ENTREZID", keytype = "ENSEMBLTRANS")
ALLt <- unique(counts_df_t[!is.na(counts_df_t$ENTREZ),]$ENTREZ)
dim(counts_df_t)  # 68722

# Selecting DE transcripts: for union(DE trs from all KO-Ctrl comps) - Ctrl-Ctrl
all_det_names <- setdiff(union(union(rownames(comp1t), rownames(comp2t)),
                              union(rownames(comp3t), rownames(comp4t))), rownames(comp5t)) # 7749
all_det_df <- as.data.frame(all_det_names)
head(all_det_df)
all_det_df$mean_lfc <- sapply(all_det_df$all_det_names, function(name) {
  ifelse(name %in% counts_df_t$tx_id, 
         counts_df_t$mean_lfc[which(counts_df_t$tx_id == name)], 
         0)
})
# Adding the entrez id 
# keytypes(org.Hs.eg.db) could be used to check which IDs are in this database
all_det_df$ENTREZ <- mapIds(org.Hs.eg.db, keys = all_det_df$all_det_names,
                           column = "ENTREZID", keytype = "ENSEMBLTRANS")
# Leaving excising entrez names + highest value among the duplicated value 
all_det_df <- all_det_df[!is.na(all_det_df$ENTREZ),]
all_det_df$mean_lfc <- as.numeric(all_det_df$mean_lfc)
all_det_df <- all_det_df %>%
  group_by(ENTREZ) %>%
  slice_max(order_by = mean_lfc, n = 1) %>%
  ungroup()
all_det_spia_list <- all_det_df$mean_lfc
names(all_det_spia_list) <- as.vector(all_det_df$ENTREZ)

# Performing SPIA analysis DE-vs-all
res_de <- spia(de=all_det_spia_list, all=ALLt, organism="hsa",
               nB=500,plots=FALSE,beta=NULL,combine="fisher",verbose=FALSE)
res_de[1:10, c(1, 11)]
write.csv(res_de, file = "../spia/spia_res_rnaseq_fullDE_transcripts.csv", row.names = TRUE)


############################ SPIA: for key genes/trs ###########################

# Saving all key proteins to the dataframe
key_genes_spia <- key_genes_gsea[, c(13, 14)]
key_genes_spia$ENTREZ <- mapIds(org.Hs.eg.db, keys = key_genes_spia$gene_name,
                                column = "ENTREZID", keytype = "SYMBOL")
key_genes_spia_list <- key_genes_spia$mean_lfc
names(key_genes_spia_list) <- as.vector(key_genes_spia$ENTREZ)

# Performing SPIA analysis key_proteins-vs-all
res <- spia(de=key_genes_spia_list, all=ALL, organism="hsa",
            nB=500,plots=FALSE,beta=NULL,combine="fisher",verbose=FALSE)
res[1:10, c(1, 11)]
write.csv(res, file = "../spia/spia_res_rnaseq_key_all.csv", row.names = TRUE)

# We can also perform separated analysis for up- and down-regulated key proteins:

# UP: data preparation
key_genes_spia_up <- subset(key_genes_gsea, gene_name %in% up_genes)[, c(13, 14)]
key_genes_spia_up$ENTREZ <- mapIds(org.Hs.eg.db, keys = key_genes_spia_up$gene_name,
                                   column = "ENTREZID", keytype = "SYMBOL")
key_genes_spia_up_list <- key_genes_spia_up$mean_lfc
names(key_genes_spia_up_list) <- as.vector(key_genes_spia_up$ENTREZ)

# UP: running SPIA
res_up <- spia(de=key_genes_spia_up_list, all=ALL, organism="hsa",
               nB=2000,plots=FALSE,beta=NULL,combine="fisher",verbose=FALSE)
res_up[1:10, c(1, 11)]
write.csv(res_up, file = "spia_res_key_up.csv", row.names = TRUE)

# DOWN: data preparation
key_genes_spia_down <- subset(key_genes_gsea, gene_name %in% down_genes)
key_genes_spia_down$ENTREZ <- mapIds(org.Hs.eg.db, keys = key_genes_spia_down$gene_name,
                                     column = "ENTREZID", keytype = "SYMBOL")
key_genes_spia_down_list <- key_genes_spia_down$mean_lfc
names(key_genes_spia_down_list) <- as.vector(key_genes_spia_down$ENTREZ)

# DOWN: running SPIA
res_down <- spia(de=key_genes_spia_down_list, all=ALL, organism="hsa",
                 nB=2000,plots=FALSE,beta=NULL,combine="fisher",verbose=FALSE)
res_down[1:10, c(1, 11)]
write.csv(res_down, file = "spia_res_key_down.csv", row.names = TRUE)


################  SPIA: for proteomics, full DE and key proteins ###############

# All proteins to compare to (background) - just preliminary filtering
df_prots <- as.data.frame(rownames(assay(cptac_imp))) # 3374
colnames(df_prots) <- c("prot_names")
df_prots$ENTREZ <- mapIds(org.Hs.eg.db, keys = df_prots$prot_names,
                             column = "ENTREZID", keytype = "SYMBOL")
ALLp <- unique(df_prots[!is.na(df_prots$ENTREZ),]$ENTREZ)
length(ALLp) # 3160


# 1) DE PROTEINS: data preparation

# Here res_AB16_imp is already filtered by logFC and p-value
prot_de_names <- union(res_AB16_imp$protein, res_AE4_imp$protein) # 1283
prot_de_df <- as.data.frame(prot_de_names)
head(prot_de_df)
# Here res_AB16_imp is NOT FILTERED!
prot_de_df$mean_lfc <- sapply(prot_de_df$prot_de_names, function(name) {
  mean(res_AB16_imp$logFC[which(res_AB16_imp$protein == name)], 
       res_AE4_imp$logFC[which(res_AE4_imp$protein == name)])
})
# Adding the entrez ID
prot_de_df$ENTREZ <- mapIds(org.Hs.eg.db, keys = prot_de_df$prot_de_names,
                           column = "ENTREZID", keytype = "SYMBOL")
# Leaving excising entrez names + highest value among the duplicated value 
prot_de_df <- prot_de_df[!is.na(prot_de_df$ENTREZ),] # 1193
prots_de_spia_list <- prot_de_df$mean_lfc
names(prots_de_spia_list) <- as.vector(prot_de_df$ENTREZ)

# 2) KEY PROTEINS from proteomics: data preparation 

key_prots_pr <- c("BRIX1", "CANX", "DICER1", "EIF4A1", "GNL3", "GAPDH", "HSP90AA1", 
                  "IDH2", "ITGA5", "KDM1A", "KIAA1429", "NOP14", "PES1", "PFDN4", 
                  "PRC1", "PRDX2", "RPL11", "RPL21", "RPS15", "RPS20", "RPS23", 
                  "TXN", "TUBB4A", "VPS52")
prot_key_df <- as.data.frame(key_prots_pr)
head(prot_key_df)
# Here res_AB16_imp is NOT FILTERED
prot_key_df$mean_lfc <- sapply(prot_key_df$key_prots_pr, function(name) {
  mean(res_AB16_imp$logFC[which(res_AB16_imp$protein == name)], 
       res_AE4_imp$logFC[which(res_AE4_imp$protein == name)])
})
# Adding the entrez ID
prot_key_df$ENTREZ <- mapIds(org.Hs.eg.db, keys = prot_key_df$key_prots_pr,
                            column = "ENTREZID", keytype = "SYMBOL")
# Leaving excising entrez names + highest value among the duplicated value 
prot_key_df <- prot_key_df[!is.na(prot_key_df$ENTREZ),] # 1193
prots_key_spia_list <- prot_key_df$mean_lfc
names(prots_key_spia_list) <- as.vector(prot_key_df$ENTREZ)


# 3) Performing SPIA analysis 

# SPIA: DE-vs-all
res_prot_de <- spia(de=prots_de_spia_list, all=ALLp, organism="hsa",
                    nB=500,plots=FALSE,beta=NULL,combine="fisher",verbose=FALSE)
res_prot_de[1:10, c(1, 11)]
write.csv(res_prot_de, file = "../spia/spia_res_prot_fullDE.csv", row.names = TRUE)

# SPIA: key-vs-all
res_prot_key <- spia(de=prots_key_spia_list, all=ALLp, organism="hsa",
                     nB=500,plots=FALSE,beta=NULL,combine="fisher",verbose=FALSE)
res_prot_key[1:10, c(1, 11)]
write.csv(res_prot_key, file = "../spia/spia_res_prot_key24.csv", row.names = TRUE)


########################## Venn diagram for SPIA terms #########################

# We can create a Venn diagram to illustrate common SPIA terms between terms that 
# we got for proteomics, and genes and transcripts from the transcriptomics data.

# PROTEINS
gene_sets <- list(DE_pathways = res_prot_de$ID, key_pathways = res_prot_key$ID)
circle_colors <- c("DE_pathways" = "darkseagreen2", "key_pathways" = "gold2")
euler_result <- euler(gene_sets, shape = "ellipse", fills = circle_colors)
plot(euler_result, quantities = TRUE, fills = circle_colors, counts = TRUE,
     legend = TRUE, main = "SPIA: DE and key proteins from proteomics")

# GENES and TRANSCRIPTS
res_de_genes <- read_csv("../spia/spia_res_rnaseq_fullDE_genes.csv")
res_de_trs <- read_csv("../spia/spia_res_rnaseq_fullDE_transcripts.csv")

gene_sets_res <- list(genes_DE = res_de_genes$ID, transcripts_DE = res_de_trs$ID,
                  key_proteins = res$ID)
circle_colors <- c("genes_DE" = "lightsalmon", "transcripts_DE" = "darkseagreen2",
                   "key_proteins" = "gold2")
euler_result <- euler(gene_sets_res, shape = "ellipse", fills = circle_colors)
plot(euler_result, quantities = TRUE, fills = circle_colors, counts = TRUE,
     legend = TRUE, main = "SPIA: DE and key proteins from transcriptomics")


###################################### GO ######################################

# Besides SPIA analysis, we can also perform Gene Onthology (GO) analysis in R
# for our results. We will also replicate this analysis using the STRING web portal. 

# GENES:
# Full list of diff. expressed genes
go <- goana(all_de_names, species="Hs", FDR = 0.1)
# selecting biological processes (BP)
top_terms <- head(filter(go[order(go$P.DE,decreasing = FALSE),], Ont=='BP'), n=100)
write.table(top_terms, "../go/top_100_GO_genes_DE.txt", row.names = FALSE,
            col.names = FALSE, quote = FALSE)

# TRANSCRIPTS:
# Full list of diff. expressed transcripts
go <- goana(all_det_names, species="Hs", FDR = 0.1)
# selecting biological processes (BP)
top_terms <- head(filter(go[order(go$P.DE,decreasing = FALSE),], Ont=='BP'), n=100)
write.table(top_terms, "../go/top_100_GO_trs_DE.txt", row.names = FALSE,
            col.names = FALSE, quote = FALSE)

# KEY GENES:
go <- goana(names(key_genes_spia_up_list), species="Hs", FDR = 0.1)
# selecting biological processes (BP)
top_terms <- head(filter(go[order(go$P.DE,decreasing = FALSE),], Ont=='BP'), n=100)
write.table(top_terms, "../go/top_100_GO_genes_key.txt", row.names = FALSE,
            col.names = FALSE, quote = FALSE)

# PROTEINS:
# Full list of diff. expressed proteins
go <- goana(names(prots_de_spia_list), species="Hs", FDR = 0.1)
# selecting biological processes (BP)
top_terms <- head(filter(go[order(go$P.DE,decreasing = FALSE),], Ont=='BP'), n=100)
write.table(top_terms, "../go/top_100_GO_prots_DE.txt", row.names = FALSE,
            col.names = FALSE, quote = FALSE)

# Only key proteins
go <- goana(names(prots_key_spia_list), species="Hs", FDR = 0.1)
# selecting biological processes (BP)
top_terms <- head(filter(go[order(go$P.DE,decreasing = FALSE),], Ont=='BP'), n=100)
write.table(top_terms, "../go/top_100_GO_prots_key.txt", row.names = FALSE,
            col.names = FALSE, quote = FALSE)


################### tRNA-related GO terms and genes ############################

# Here we are extracting genes related to tRNA-processing GO terms. These terms
# were extracted from the Gene Onthology portal and were saved to 'tRNA_all_go_terms.tsv'.
# There are 185 such GO terms in the file.

# Uploading my GO list
go_list <- read.csv("/home/anna/Downloads/OE33_project/trna_prots/tRNA_all_go_terms.csv",
                    sep='\t', header = FALSE)
colnames(go_list) <- c("GO_id", "GO_name")

# In this loop we are iterating over the list of tRNA-related terms and selecting
# genes from the org.Hs.eg.db database, assosiated with these GO terms.

genes_list <- c()
for(go in go_list$GO_id) {
  # Use tryCatch to handle potential errors
  tryCatch({
    genes <- AnnotationDbi::select(org.Hs.eg.db, go, c("GENENAME","SYMBOL"), c("GO"))
    genes_list <- append(genes_list, genes$SYMBOL)
  }, error = function(e) {
    cat(paste("Can't find any human genes for",go,":(\n")) # conditionMessage(e)
  })
}

genes_list <- unique(genes_list) # 256 unique genes overall

# Additionally, GO:0006399 tRNA metabolic process is mentioned on the Gene Onthology 
# portal as an umbrella GO term. So we extracted genes related to this term from the portal.
# There are 199 overal, 181 of which are in genes_list.
# We will keep the union of two lists as a final list of tRNA-related genes.

go_fr_one_term <- read.table("/home/anna/Downloads/OE33_project/trna_prots/tRNA_go_one_term.txt",
                             header = FALSE)$V1
genes_list <- union(genes_list, go_fr_one_term)
write.table(genes_list, file = "/home/anna/Downloads/OE33_project/trna_prots/tRNA_all_GO_genes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# Moreover, the publication mining showed there are seleral proteins potentially
# involved in other tRNA chemical modifications that increase tRNA stability. 
# There are 30 of such genes, and 21 of them are in the list mentioned above. 

goi_pub_list <- read.table("/home/anna/Downloads/OE33_project/trna_prots/poi.txt", header = FALSE)$V1
length(intersect(goi_pub_list, genes_list)) # 21

# We can check then how many of these proteins are in our 'general intersection'
# derived from transcriptomics or in the proteomics dataset. 
filter(inters_AB16_AE4_imp_vals, prnames %in% intersect(genes_list, inters_AB16_AE4_imp_names))
filter(global_int_t_counts, gene_name %in% intersect(genes_list, global_int_t_genenames))

########################## Saving intersection data for GSEA ###################

# One more type of functional analysis we perform is the gene set enrichment 
# analysis (GSEA). In order to run it in the application, we have to prepare
# data in a certain format. 

###### for the full dataset:

norm_counts <- counts(dds, normalized = T)
norm_counts <- as.data.frame(norm_counts)
norm_counts$NAME <- rownames(norm_counts)
norm_counts$Description <- rownames(norm_counts)
new_order <- c('NAME', 'Description', setdiff(names(norm_counts), c('NAME', 'Description')))
norm_counts <- norm_counts[, new_order]

fid <- "gsea_data_full.gct"
writeLines(c("#1.2", paste(nrow(norm_counts),
                           ncol(norm_counts) - 2, collapse="\t")), fid, sep="\n")
write.table(norm_counts, file=fid, quote=FALSE, row.names=FALSE,
            col.names=TRUE, sep="\t", append = TRUE)

####### for key proteins: 

key_genes_gsea <- subset(global_int_g_counts, gene_name %in% key_prots)

# up-regulated
key_genes_gsea_up <- subset(key_genes_gsea, gene_name %in% up_genes)
key_genes_gsea_up$NAME <- key_genes_gsea_up$gene_name
key_genes_gsea_up$Description <- key_genes_gsea_up$gene_name
key_genes_gsea_up$gene_name <- NULL
key_genes_gsea_up$mean_lfc <- NULL
new_order <- c('NAME', 'Description', setdiff(names(key_genes_gsea_up), c('NAME', 'Description')))
key_genes_gsea_up <- key_genes_gsea_up[, new_order]

fid <- "gsea_data_key_up.gct"
writeLines(c("#1.2", paste(nrow(key_genes_gsea_up),
                           ncol(key_genes_gsea_up) - 2, collapse="\t")), fid, sep="\n")
write.table(key_genes_gsea_up, file=fid, quote=FALSE, row.names=FALSE,
            col.names=TRUE, sep="\t", append = TRUE)
write.table(key_genes_gsea_up$NAME, file='up_key_names.txt', quote=FALSE, row.names=FALSE,
            col.names=F)

# down-regulated
key_genes_gsea_down <- subset(key_genes_gsea, gene_name %in% down_genes)
key_genes_gsea_down$NAME <- key_genes_gsea_down$gene_name
key_genes_gsea_down$Description <- key_genes_gsea_down$gene_name
key_genes_gsea_down$gene_name <- NULL
key_genes_gsea_down$mean_lfc <- NULL
new_order <- c('NAME', 'Description', setdiff(names(key_genes_gsea_down), c('NAME', 'Description')))
key_genes_gsea_down <- key_genes_gsea_down[, new_order]
fid <- "gsea_data_key_down.gct"
writeLines(c("#1.2", paste(nrow(key_genes_gsea_down),
                           ncol(key_genes_gsea_down) - 2, collapse="\t")),
           fid, sep="\n")
write.table(key_genes_gsea_down, file=fid, quote=FALSE, row.names=FALSE,
            col.names=TRUE, sep="\t", append = TRUE)
