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
library(eulerr)

setwd("~/Downloads/OE33_project")

############### Intersection of EdgeR-DESeq2 results (names) ###################

# Gene/trs names from edgeR are in: filtered_contrasts_g_names, filtered_contrasts_t_names
# Gene/trs names from DESeq2 are in: deseq_filtered_g_names, deseq_filtered_t_names
# Names of comparisons are in celllines_pairs

# Plotting
circle_colors <- c("DESeq2" = "skyblue", "edgeR" = "lightcoral")

# GENES:
dev.off()
for(i in 1:5){
  gene_sets <- list(DESeq2 = deseq_filtered_g_names[[i]],
                    edgeR = filtered_contrasts_g_names[[i]])
  euler_result <- euler(gene_sets, shape = "ellipse", fills = circle_colors)
  show(plot(euler_result, quantities = TRUE, fills = circle_colors, counts = TRUE,
       legend = TRUE, main = paste("Genes:", gsub("_", " ", celllines_pairs[[i]]))))
}
# TRANSCRIPTS:
dev.off()
for(i in 1:5){
  gene_sets <- list(DESeq2 = deseq_filtered_t_names[[i]],
                    edgeR = filtered_contrasts_t_names[[i]])
  euler_result <- euler(gene_sets, shape = "ellipse", fills = circle_colors)
  show(plot(euler_result, quantities = TRUE, fills = circle_colors, counts = TRUE,
       legend = TRUE, main = paste("Trs:", gsub("_", " ", celllines_pairs[[i]]))))
}

# Finding intersections between gene/transcrsipts names
gene_intersections <- c('ae4_ab1_intg', 'ab16_ab1_intg', 'ae4_ae21_intg',
                        'ab16_ae21_intg', 'ab1_ae21_intg')
transcript_intersections <- c('ae4_ab1_intt', 'ab16_ab1_intt', 'ae4_ae21_intt',
                              'ab16_ae21_intt', 'ab1_ae21_intt')

ints_genes_list <- list()
ints_trs_list <- list()

for (i in 1:5) {
  ints_genes_list[[gene_intersections[i]]] <- intersect(deseq_filtered_g_names[[i]], filtered_contrasts_g_names[[i]])
  ints_trs_list[[transcript_intersections[i]]] <- intersect(deseq_filtered_t_names[[i]], filtered_contrasts_t_names[[i]])
}

# Saving names
# GENES:
output_dir <- './edger_deseq_inters_all_comp'
type_dir <- 'genes'
for(df_name in names(ints_genes_list)){
  write.table(ints_genes_list[df_name], file = file.path(output_dir, type_dir, paste0("genes_",df_name,"_names.txt")),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}
# TRANSCRIPTS:
type_dir <- 'transcripts'
for(df_name in names(ints_trs_list)){
  write.table(ints_trs_list[df_name], file = file.path(output_dir, type_dir, paste0("trs_",df_name,"_names.txt")),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

######### Saving the EdgeR-DESeq2 intersections (with EdgeR stats) #############

# We will use these sets for the primary network analysis, here we are keeping
# genes, that show diff. expression between two controls (ctrl-ctrl).

ints_genes_full_list <- list()
ints_trs_full_list <- list()

for (i in 1:5) {
  compname_g <- gene_intersections[i]
  compname_t <- transcript_intersections[i]
  ints_genes_full_list[[compname_g]] <- subset(filtered_contrasts_g[[i]][, c(1, 4, 6)],
                                              gene_name %in% ints_genes_list[[compname_g]])
  ints_trs_full_list[[compname_t]] <- subset(filtered_contrasts_t[[i]][, c(1, 4, 6)],
                                                  rownames(filtered_contrasts_t[[i]][, c(1, 4, 6)]) %in% ints_trs_list[[compname_t]])
}

# GENES:
type_dir <- 'genes'
for(df_name in names(ints_genes_full_list)){
  write.csv(ints_genes_full_list[df_name], file = file.path(output_dir, type_dir,paste0("genes_",df_name,"_full.txt")),
              row.names = TRUE)
}
# TRANSCRIPTS:
type_dir <- 'transcripts'
for(df_name in names(ints_trs_full_list)){
  write.csv(ints_trs_full_list[df_name], file = file.path(output_dir, type_dir, paste0("trs_",df_name,"_full.txt")),
            row.names = TRUE)
  write.table(ints_trs_full_list[df_name]$gene_name, file = file.path(output_dir, type_dir, paste0("trs_",df_name,"_gene_names.txt")),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}


####### 2nd level subtracting ctrl-ctrl, intersections, venn diagrams ##########

# Here we are performing the next step of intersections, by intersecting samples
# with same background (e.g. what is common between KO1-vs-Ctrl1 and KO2-vs-Ctrl1), 
# and with same experimental samples (e.g. what is common between KO1-vs-Ctrl1 and KO1-vs-Ctrl2).
# By doing this we can extract the most consistently diff. expressed sets of genes/transcripts. 

# Function for plotting venn diagram for sets of genes/trs names
plot_inters_2lvl <- function(gene_sets, gene_or_trs){
  names <- names(gene_sets)
  circle_colors <- setNames(c('lightgreen', 'yellow', 'coral'), names)
  euler_result <- euler(gene_sets, shape = "ellipse", fills = circle_colors)
  show(plot(euler_result, quantities = TRUE, fills = circle_colors, counts = TRUE,
       legend = TRUE, main = paste(gene_or_trs,names[1],"vs",names[2],"\n+ AB1-AE21")))
}

# GENES:
gene_set1 <- list('AE4-AB1' = ints_genes_list[['ae4_ab1_intg']],
                   'AB16-AB1' = ints_genes_list[['ab16_ab1_intg']],
                   'AB1-AE21' = ints_genes_list[['ab1_ae21_intg']])
gene_set2 <- list('AE4-AE21' = ints_genes_list[['ae4_ae21_intg']],
                   'AB16-AE21' = ints_genes_list[['ab16_ae21_intg']],
                   'AB1-AE21' = ints_genes_list[['ab1_ae21_intg']])
gene_set3 <- list('AE4-AB1' = ints_genes_list[['ae4_ab1_intg']],
                   'AE4-AE21' = ints_genes_list[['ae4_ae21_intg']],
                   'AB1-AE21' = ints_genes_list[['ab1_ae21_intg']])
gene_set4 <- list('AB16-AB1' = ints_genes_list[['ab16_ab1_intg']],
                   'AB16-AE21' = ints_genes_list[['ab16_ae21_intg']],
                   'AB1-AE21' = ints_genes_list[['ab1_ae21_intg']])

for(set in list(gene_set1,gene_set2,gene_set3,gene_set4)){
  plot_inters_2lvl(set, "Genes:")
}

# TRANSCRIPTS:
tr_set1 <- list('AE4-AB1' = ints_trs_list[['ae4_ab1_intt']],
                    'AB16-AB1' = ints_trs_list[['ab16_ab1_intt']],
                    'AB1-AE21' = ints_trs_list[['ab1_ae21_intt']])
tr_set2 <- list('AE4-AE21' = ints_trs_list[['ae4_ae21_intt']],
                    'AB16-AE21' = ints_trs_list[['ab16_ae21_intt']],
                    'AB1-AE21' = ints_trs_list[['ab1_ae21_intt']])
tr_set3 <- list('AE4-AB1' = ints_trs_list[['ae4_ab1_intt']],
                    'AE4-AE21' = ints_trs_list[['ae4_ae21_intt']],
                    'AB1-AE21' = ints_trs_list[['ab1_ae21_intt']])
tr_set4 <- list('AB16-AB1' = ints_trs_list[['ab16_ab1_intt']],
                    'AB16-AE21' = ints_trs_list[['ab16_ae21_intt']],
                    'AB1-AE21' = ints_trs_list[['ab1_ae21_intt']])

for(set in list(tr_set1,tr_set2,tr_set3,tr_set4)){
  plot_inters_2lvl(set, "Trs:")
}

# Saving these sets
same_back_g_ab1 <- setdiff(intersect(gene_set1[['AE4-AB1']], gene_set1[['AB16-AB1']]),
                           gene_set1[['AB1-AE21']])
same_back_g_ae21 <- setdiff(intersect(gene_set2[['AE4-AE21']], gene_set2[['AB16-AE21']]),
                            gene_set2[['AB1-AE21']])
same_ko_g_ae4 <- setdiff(intersect(gene_set3[['AE4-AB1']], gene_set3[['AE4-AE21']]),
                         gene_set3[['AB1-AE21']])
same_ko_g_ab16 <- setdiff(intersect(gene_set4[['AB16-AB1']], gene_set4[['AB16-AE21']]),
                          gene_set4[['AB1-AE21']])

same_back_t_ab1 <- setdiff(intersect(tr_set1[['AE4-AB1']], tr_set1[['AB16-AB1']]),
                           tr_set1[['AB1-AE21']])
same_back_t_ae21 <- setdiff(intersect(tr_set2[['AE4-AE21']], tr_set2[['AB16-AE21']]),
                            tr_set2[['AB1-AE21']])
same_ko_t_ae4 <- setdiff(intersect(tr_set3[['AE4-AB1']], tr_set3[['AE4-AE21']]),
                         tr_set3[['AB1-AE21']])
same_ko_t_ab16 <- setdiff(intersect(tr_set4[['AB16-AB1']], tr_set4[['AB16-AE21']]),
                          tr_set4[['AB1-AE21']])


##### 3nd level intersections = "global intersection", saving with stats #######

# Plotting the global intersection - intersection of 2nd level intersections that
# we created in the previous section. 

# GENES:
gene_sets <- list('KOs-AB1' = same_back_g_ab1, 'KOs-AE21' = same_back_g_ae21,
                  'AE4-Ctrls' =  same_ko_g_ae4, 'AB16-Ctrls' = same_ko_g_ab16)
circle_colors <- c('KOs-AB1' = 'green', 'KOs-AE21' = 'yellow',
                   'AE4-Ctrls' =  'coral', 'AB16-Ctrls' = 'blue')
euler_result <- euler(gene_sets, shape = "ellipse", fills = circle_colors)
plot(euler_result, quantities = TRUE, fills = circle_colors, counts = TRUE,
     legend = TRUE, main = "Genes: all intercestions w/o Ctlr-Ctrl")

# TRANSCRIPTS:
gene_sets <- list('KOs-AB1' = same_back_t_ab1, 'KOs-AE21' = same_back_t_ae21,
                  'AE4-Ctrls' =  same_ko_t_ae4, 'AB16-Ctrls' = same_ko_t_ab16)
circle_colors <- c('KOs-AB1' = 'green', 'KOs-AE21' = 'yellow',
                   'AE4-Ctrls' =  'coral', 'AB16-Ctrls' = 'blue')
euler_result <- euler(gene_sets, shape = "ellipse", fills = circle_colors)
plot(euler_result, quantities = TRUE, fills = circle_colors, counts = TRUE,
     legend = TRUE, main = "Transcripts: all intercestions w/o Ctlr-Ctrl")

# Saving the 'global intersection' to dataframes
global_int_g <- intersect(intersect(same_back_g_ab1, same_back_g_ae21),
                         intersect(same_ko_g_ae4, same_ko_g_ab16))
global_int_t <- intersect(intersect(same_back_t_ab1, same_back_t_ae21),
                         intersect(same_ko_t_ae4, same_ko_t_ab16))

# Now we need to filter 'global intersection' results so all genes from all comparisons 
# are diff. expressed in ONE DIRECTION, e.g. in ALL comparisons KO-CTRL we have either
# logFC < 0.5 or logFC > 0.5.

# Supporting function to check consistency in direction over edgeR datasets
check_dir <- function(name, type){
  if(type == 'trs'){
    comps_list <- list('c1' = comp1t, 'c2' = comp2t, 'c3' = comp3t, 'c4' = comp4t)
    if(all(sapply(comps_list, function(comp) all(subset(comp, rownames(comp) == name)$logFC > 0)))){
      return('all_up')
    } else if(all(sapply(comps_list, function(comp) all(subset(comp, rownames(comp) == name)$logFC < 0)))){
      return('all_down')
    }
  } else {
    comps_list <- list('c1' = comp1g, 'c2' = comp2g, 'c3' = comp3g, 'c4' = comp4g)
    if(all(sapply(comps_list, function(comp) all(subset(comp, gene_name == name)$logFC > 0)))){
      return('all_up')
    } else if(all(sapply(comps_list, function(comp) all(subset(comp, gene_name == name)$logFC < 0)))){
      return('all_down')
    }
  }
  return('none')
}

# Function that takes names of genes/trs and returns list of names of consistently
# up- and down-regulated genes/trs
select_one_direction <- function(inters_names){
  if(substr(inters_names[1], 4, 4) == "T"){
    type <- 'trs'
  } else {
    type <- 'genes'
  }
  down_names <- c()
  up_names <- c()
  for(name in inters_names){
    if(check_dir(name, type) == "all_up"){
      up_names <- append(up_names, name)
    } else {
      if(check_dir(name, type) == "all_down"){
        down_names <- append(down_names, name)
      }
    }
  }
  return(list('down' = down_names, 'up' = up_names))
}

# Function to calculate mean logFC values for given genes over several comparisons
calc_mean_lfc <- function(data, type){
  new_data <- data
  new_data$ids <- rownames(data)
  if(type == 'trs'){
    comp1t$ids <- rownames(comp1t)
    new_data$tx_name <- NULL
    # Adding gene_names column
    new_data <- merge(new_data, comp1t[,6:7], by.x = 'ids', by.y = 'ids')
    # Adding mean logFC column for each id=row
    new_data$mean_lfc <- sapply(new_data$ids, function(name) {
      mean(sapply(filtered_contrasts_t[1:4], function(df) df[name, "logFC"]))
    })
  } else {
    filtered_contrasts_g[[1]]$ids <- rownames(filtered_contrasts_g[[1]])
    # Adding mean logFC column for each id=row
    new_data$mean_lfc <- sapply(new_data$ids, function(name) {
      mean(sapply(filtered_contrasts_g[1:4], function(df) df[name, "logFC"]))
    })
  }
  rownames(new_data) <- new_data$ids
  new_data$ids <- NULL
  return(new_data)
}

# Marking genes/trs having all logFC in the same direction
down_genes <- select_one_direction(global_int_g)[['down']]
up_genes <- select_one_direction(global_int_g)[['up']]
down_trs <- select_one_direction(global_int_t)[['down']]
up_trs <- select_one_direction(global_int_t)[['up']]

# Merging up- and down-regulated genes/trs
global_int_g_filtered <- append(down_genes, up_genes)
global_int_t_filtered <- append(down_trs, up_trs)
length(global_int_g_filtered)

global_int_t_genenames <- subset(ae4_ab1_intt_full, rownames(ae4_ab1_intt_full) %in% global_int_t_filtered)$gene_name
tr_gene_inters <- intersect(global_int_t_genenames, global_int_g_filtered)
length(tr_gene_inters)

# Extracting raw counts for these genes/trs
counts_df_g <- as.data.frame(fit_genes$counts)
counts_df_g$gene_name <- fit_genes$genes$gene_name
global_int_g_counts <- subset(counts_df_g, gene_name %in% global_int_g_filtered)
global_int_g_counts <- calc_mean_lfc(global_int_g_counts, 'gene')

counts_df_t <- as.data.frame(fit_tr$counts)
counts_df_t$tx_name <- fit_tr$genes$tx
global_int_t_counts <- subset(counts_df_t, tx_name %in% global_int_t_filtered)
global_int_t_counts <- calc_mean_lfc(global_int_t_counts, 'trs')

# Saving 'global intersection'

type_dir <- 'global_intersection'

# GENES: names and counts
write.csv(global_int_g_counts, file = file.path(output_dir, type_dir, "Global_int_genes_counts.csv"), row.names = TRUE)
write.table(global_int_g_filtered, file = file.path(output_dir, type_dir, "Global_int_genes_names.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# TRANSCRIPTS: names and counts
write.csv(global_int_t_counts, file = file.path(output_dir, type_dir, "Global_int_transcripts_counts.csv"), row.names = TRUE)
write.table(global_int_t_filtered, file = file.path(output_dir, type_dir, "Global_int_transcripts_names.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(global_int_t_genenames, file = file.path(output_dir, type_dir, "Global_int_transcripts_genenames.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
# GENES-TRS INTERSECTION: gene names
write.table(tr_gene_inters, file = file.path(output_dir, type_dir, "Global_int_trs_genes_inters.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)


####################### Heatmaps with 'global intersection' ####################

######### Transcripts, all

global_int_t_counts_hm <- global_int_t_counts[, 1:12]
head(global_int_t_counts_hm)
names(global_int_t_counts_hm) <- c("Ctrl_AB1_1", "Ctrl_AB1_2", "Ctrl_AB1_3", 
                                   "Ctrl_AE21_1", "Ctrl_AE21_2", "Ctrl_AE21_3",
                                   "KO_AB16_1", 'KO_AB16_2', 'KO_AB16_3', 
                                   'KO_AE4_1', 'KO_AE4_2', 'KO_AE4_3')
# For col clusters
samples <- data.frame(sample = as.character(colnames(global_int_t_counts_hm)), Sample_type="METTL1_KO") %>%
  column_to_rownames("sample")
samples$Sample_type <- c(rep("Ctrl", 6), rep("METTL1_KO", 6))
# For row clusters
up_down_genes <- data.frame(row = as.character(rownames(global_int_t_counts_hm)), Regulation = "Up")
rownames(up_down_genes) <- up_down_genes$row
up_down_genes$Regulation <- ifelse(rownames(up_down_genes) %in% down_trs, "Down", "Up")

o <- order(up_down_genes$Regulation)
up_down_genes <- up_down_genes[o, ]
global_int_t_counts_hm <- global_int_t_counts_hm[o, ]
up_down_genes$row <- NULL

hm_gg <- as.ggplot(pheatmap(as.matrix(log(global_int_t_counts_hm+1)), scale="row",
                            cluster_rows = F,
                            fontsize_col = 8,
                            #fontsize_row = 3,
                            angle_col = 90,
                            cellheight = 2,
                            treeheight_col = 3,
                            treeheight_row = 10,
                            annotation_row = up_down_genes,
                            annotation_col= samples,
                            annotation_colors =list(
                              Sample_type=c(METTL1_KO="grey",Ctrl="black"),
                              Regulation=c(Up="chocolate1", Down="cornflowerblue")),
                            color=colorRampPalette(c("navy","white", "red"))(50), 
                            show_rownames = F, annotation_names_col = F, annotation_names_row = F))


######### Genes, top-30 
# All the 512 genes are not fitting, so we are selecting ones with abs(logFC) > 1

global_int_g_counts_hm <- global_int_g_counts
global_int_g_counts_hm_top30 <- head(
  global_int_g_counts_hm[order(abs(global_int_g_counts_hm$mean_lfc), decreasing = T),],
  n=30)
names(global_int_g_counts_hm_top30) <- c("Ctrl_AB1_1", "Ctrl_AB1_2", "Ctrl_AB1_3", 
                                         "Ctrl_AE21_1", "Ctrl_AE21_2", "Ctrl_AE21_3",
                                         "KO_AB16_1", 'KO_AB16_2', 'KO_AB16_3', 
                                         'KO_AE4_1', 'KO_AE4_2', 'KO_AE4_3', 'gene_name', 
                                         'mean_lfc')

# For row clusters
up_down_genes <- data.frame(row = as.character(rownames(global_int_g_counts_hm_top30)),
                            gene_name = as.character(global_int_g_counts_hm_top30$gene_name))
rownames(up_down_genes) <- up_down_genes$gene_name
up_down_genes$Regulation <- ifelse(up_down_genes$gene_name %in% down_genes, "Down", "Up")

o <- order(up_down_genes$Regulation)
up_down_genes <- up_down_genes[o, ]
rownames(global_int_g_counts_hm_top30) <- global_int_g_counts_hm_top30$gene_name
global_int_g_counts_hm_top30 <- global_int_g_counts_hm_top30[o, 1:12]
up_down_genes$row <- NULL
up_down_genes$gene_name <- NULL

hm_gg <- as.ggplot(pheatmap(as.matrix(log(global_int_g_counts_hm_top30+1)), scale="row",
                            cluster_rows = T,
                            fontsize_col = 8,
                            #fontsize_row = 3,
                            angle_col = 90,
                            cellheight = 10,
                            treeheight_col = 0,
                            treeheight_row = 0,
                            annotation_row = up_down_genes,
                            annotation_col= samples,
                            annotation_colors =list(
                              Sample_type=c(METTL1_KO="grey",Ctrl="black"),
                              Regulation=c(Up="chocolate1", Down="cornflowerblue")),
                            color=colorRampPalette(c("navy","white", "red"))(50), 
                            show_rownames = T, annotation_names_col = F, annotation_names_row = F))


######### Heatmap for top-30 of TRANSCRIPTS

global_int_t_counts_hm <- global_int_t_counts
global_int_t_counts_hm$tx_name <- NULL
global_int_t_counts_hm_top30 <- head(
  global_int_t_counts_hm[order(abs(global_int_t_counts_hm$mean_lfc), decreasing = T),],
  n=30)
names(global_int_t_counts_hm_top30) <- c("Ctrl_AB1_1", "Ctrl_AB1_2", "Ctrl_AB1_3", 
                                         "Ctrl_AE21_1", "Ctrl_AE21_2", "Ctrl_AE21_3",
                                         "KO_AB16_1", 'KO_AB16_2', 'KO_AB16_3', 
                                         'KO_AE4_1', 'KO_AE4_2', 'KO_AE4_3', 'gene_name', 
                                         'mean_lfc')

# for row clusters
up_down_genes <- data.frame(row = as.character(rownames(global_int_t_counts_hm_top30)),
                            gene_name = as.character(global_int_t_counts_hm_top30$gene_name))
rownames(up_down_genes) <- up_down_genes$gene_name
up_down_genes$Regulation <- ifelse(up_down_genes$row %in% down_trs, "Down", "Up")

o <- order(up_down_genes$Regulation)
up_down_genes <- up_down_genes[o, ]
rownames(global_int_t_counts_hm_top30) <- global_int_t_counts_hm_top30$gene_name
global_int_t_counts_hm_top30 <- global_int_t_counts_hm_top30[o, 1:12]
up_down_genes$row <- NULL
up_down_genes$gene_name <- NULL

hm_gg <- as.ggplot(pheatmap(as.matrix(log(global_int_t_counts_hm_top30+1)), scale="row",
                            cluster_rows = T,
                            fontsize_col = 8,
                            #fontsize_row = 3,
                            angle_col = 90,
                            cellheight = 10,
                            treeheight_col = 0,
                            treeheight_row = 0,
                            annotation_row = up_down_genes,
                            annotation_col= samples,
                            annotation_colors =list(
                              Sample_type=c(METTL1_KO="grey",Ctrl="black"),
                              Regulation=c(Up="chocolate1", Down="cornflowerblue")),
                            color=colorRampPalette(c("navy","white", "red"))(50), 
                            show_rownames = T, annotation_names_col = F, annotation_names_row = F))


######### Key GENES heatmap

# Key genes were obtained as described in methodology and then this list was
# used here to make a heatplot. 

key_prots <- read.table("./edger_deseq_inters_all_comp/global_intersection/key_prots_genenames.txt")$V1
global_int_key_genes <- subset(global_int_g_counts, gene_name %in% key_prots)

names(global_int_key_genes) <- c("Ctrl_AB1_1", "Ctrl_AB1_2", "Ctrl_AB1_3", 
                                 "Ctrl_AE21_1", "Ctrl_AE21_2", "Ctrl_AE21_3",
                                 "KO_AB16_1", 'KO_AB16_2', 'KO_AB16_3', 
                                 'KO_AE4_1', 'KO_AE4_2', 'KO_AE4_3', 'gene_name', 
                                 'mean_lfc')

# For row clusters
up_down_genes <- data.frame(row = as.character(rownames(global_int_key_genes)),
                            gene_name = as.character(global_int_key_genes$gene_name))
rownames(up_down_genes) <- up_down_genes$gene_name
up_down_genes$Regulation <- ifelse(up_down_genes$gene_name %in% down_genes, "Down", "Up")

o <- order(up_down_genes$Regulation)
up_down_genes <- up_down_genes[o, ]
rownames(global_int_key_genes) <- global_int_key_genes$gene_name
global_int_key_genes <- global_int_key_genes[o, 1:12]
up_down_genes$row <- NULL
up_down_genes$gene_name <- NULL

hm_gg <- as.ggplot(pheatmap(as.matrix(log(global_int_key_genes+1)), scale="none",
                            cluster_rows = F,
                            fontsize_col = 8,
                            fontsize_row = 6,
                            angle_col = 90,
                            cellheight = 7,
                            treeheight_col = 0,
                            treeheight_row = 0,
                            annotation_row = up_down_genes,
                            annotation_col= samples,
                            annotation_colors =list(
                              Sample_type=c(METTL1_KO="grey",Ctrl="black"),
                              Regulation=c(Up="chocolate1", Down="cornflowerblue")),
                            color=colorRampPalette(c("navy","white", "red"))(50), 
                            show_rownames = T, annotation_names_col = F, annotation_names_row = F))


######### Key TRANSCRIPTS heatmap, FILTERED

global_int_key_trs <- subset(counts_transcripts, gene_name %in% key_prots)
# If we select all transcripts from unfiltered data, that are coding for key 
# proteins, there will be 812 transcripts. All these genes are not fitting into
# the heatmap properly, so we will later make another plot to explore the
# distribution of counts among transcripts for each gene.

counts_df_t <- as.data.frame(fit_tr$counts)
counts_df_t$tx_name <- fit_tr$genes$tx
global_int_key_trs <- calc_mean_lfc(counts_df_t, 'trs')
global_int_key_trs <- subset(global_int_t_counts, gene_name %in% key_prots)
global_int_key_trs$tx_name <- NULL

names(global_int_key_trs) <- c("Ctrl_AB1_1", "Ctrl_AB1_2", "Ctrl_AB1_3", 
                               "Ctrl_AE21_1", "Ctrl_AE21_2", "Ctrl_AE21_3",
                               "KO_AB16_1", 'KO_AB16_2', 'KO_AB16_3', 
                               'KO_AE4_1', 'KO_AE4_2', 'KO_AE4_3', 'gene_name', 
                               'mean_lfc')

# For row clusters
up_down_genes <- data.frame(row = as.character(rownames(global_int_key_trs)),
                            gene_name = as.character(global_int_key_trs$gene_name))
rownames(up_down_genes) <- up_down_genes$gene_name
up_down_genes$Regulation <- ifelse(up_down_genes$row %in% down_trs, "Down", "Up")

o <- order(up_down_genes$Regulation)
up_down_genes <- up_down_genes[o, ]
rownames(global_int_key_trs) <- global_int_key_trs$gene_name
global_int_key_trs <- global_int_key_trs[o, 1:12]
up_down_genes$row <- NULL
up_down_genes$gene_name <- NULL

hm_gg <- as.ggplot(pheatmap(as.matrix(log(global_int_key_trs+1)), scale="none",
                            cluster_rows = F,
                            cluster_cols = F,
                            fontsize_col = 8,
                            fontsize_row = 8,
                            angle_col = 90,
                            cellheight = 10,
                            treeheight_col = 0,
                            treeheight_row = 0,
                            annotation_row = up_down_genes,
                            annotation_col= samples,
                            annotation_colors =list(
                              Sample_type=c(METTL1_KO="grey",Ctrl="black"),
                              Regulation=c(Up="chocolate1", Down="cornflowerblue")),
                            color=colorRampPalette(c("navy","white", "red"))(50), 
                            show_rownames = T, annotation_names_col = F, annotation_names_row = F))

