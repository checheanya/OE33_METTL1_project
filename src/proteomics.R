############################### Loading packages ###############################

library(limma)
library(dplyr)
library(reshape2)
library(stringr)
library(tibble)
library(tidyverse)
library(tidyr)
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

library(imputeLCMD)
library(MSnbase)
library(SummarizedExperiment)
library(QFeatures)

# setting directory
setwd("~/Downloads/OE33_project/Proteomics_20230721d_JR_HS_METTL1_KO")

############################# Loading pre-processed data #######################

# Description of the raw files from which the peaks data was extracted
testfile <- "./quantms_out/20230802a_JR_HS_METTL1_KO_experimental_design_mzML_openms_design_openms.mzTab"
prot <- readMzTabData(testfile, "PRT")
experimentData(prot)

# Loading intensity data for each protein per sample
intens <- read.csv("intensity_data.csv")
dim(intens)   # 38 samples 3932 proteins
length(unique(intens$accession)) # 3932 unique protein ids
length(unique(intens$symbol))  # 3869 unique protein symbols
intens %>% filter(rowSums(intens[, 3:38]) == 0) # there are no rows with all zeros

# Loading ratios after external statistical computation using contrast matrix
# (long format)
ratio <- read.csv("ratio_data.csv")
length(unique(ratio$Protein)) # 3286 unique proteins
length(unique(ratio$Symbol))  # 3282 unique protein symbols

# Statistical test results for comparisons, reshaped to the wide format
raw_stats <- read.csv("./quantms_out/20230802a_JR_HS_METTL1_KO_experimental_design_mzML_openms_design_msstats_in_comparisons.csv", 
                      sep = "\t")
length(unique(raw_stats$Protein)) # 3583 unique protein ids

# Loading a raw table with peptides and proteins annotaion
stats_not_comp <- read.csv("./quantms_out/20230802a_JR_HS_METTL1_KO_experimental_design_mzML_openms_design_msstats_in.csv")
dim(stats_not_comp) # 544857 12
length(unique(stats_not_comp$ProteinName)) # 3932 unique proteins

######################### PCA and samples clustering ###########################

non0mins <- c()
for(col in colnames(intens[,3:38])){
  non0mins <- append(non0mins, log(unique(sort(intens[,col]))[2]))
}
non0mins <- as.data.frame(non0mins)
rownames(non0mins) <- colnames(intens[,3:38])
colnames(non0mins) <- c('value')

# plotting the raw data - MANY ZEROS
ggplot(gather(intens[, 3:38]), aes(value)) + 
  geom_histogram(bins = 10) + 
  facet_wrap(~key, scales = 'free_x')

# normalizing the intensity counts - log transform
long_f_log_wt <- gather(log(as.data.frame(intens[,3:11] +1)))
long_f_log_ko_ab16 <- gather(log(as.data.frame(intens[,12:20] +1)))
long_f_log_ko_ae4 <- gather(log(as.data.frame(intens[,21:29] +1)))
long_f_log_nc2_wt <- gather(log(as.data.frame(intens[,30:38] +1)))
head(intens[,30:38])
dim(long_f_log_wt)

create_vec <- function(vec, n){
  cout <- c()
  for(val in vec){
    cout <- append(cout, rep(val, n))
  }
  return(cout)
}

long_f_log_wt$min <- create_vec(non0mins[1:9,], 3932)
long_f_log_ko_ab16$min <- create_vec(non0mins[10:18,], 3932)
long_f_log_ko_ae4$min <- create_vec(non0mins[19:27,], 3932)
long_f_log_nc2_wt$min <- create_vec(non0mins[28:36,], 3932)

# plotting the nor-ed counts
ggplot(long_f_log_wt, aes(value)) + 
  geom_histogram(bins = 30) + 
  facet_wrap(~key, scales = 'free_x') +
  geom_vline(aes(xintercept = min, color = 'red'), linetype = 'dashed') +
  xlab("log intensity") + ylab("peptide count")+
  theme(legend.position = "none")

ggplot(long_f_log_ko_ab16, aes(value)) + 
  geom_histogram(bins = 30) + 
  facet_wrap(~key, scales = 'free_x')+
  xlab("log intensity") + ylab("peptide count")+
  geom_vline(aes(xintercept = min, color = 'red'), linetype = 'dashed')+
  theme(legend.position = "none")

ggplot(long_f_log_ko_ae4, aes(value)) + 
  geom_histogram(bins = 30) + 
  facet_wrap(~key, scales = 'free_x')+
  xlab("log intensity") + ylab("peptide count")+
  geom_vline(aes(xintercept = min, color = 'red'), linetype = 'dashed')+
  theme(legend.position = "none")

ggplot(long_f_log_nc2_wt, aes(value)) + 
  geom_histogram(bins = 30) + 
  facet_wrap(~key, scales = 'free_x')+
  xlab("log intensity") + ylab("peptide count")+
  geom_vline(aes(xintercept = min, color = 'red'), linetype = 'dashed')+
  theme(legend.position = "none")

# deleting rows with all zeros
all_log <- log(intens[,3:38]+1)
# the last one - x36 is an outlier
all_log <- log(intens[,3:37]+1)
# removing rows with all zeros
all_log <- subset(all_log, rowSums(all_log != 0) > 0)

head(all_log)
pc <- prcomp(t(as.matrix(all_log)), center = TRUE, scale = TRUE)
attributes(pc)
pc$scale
summary(pc)

# Access PCA results
pca_data <- as.data.frame(pc$x)
tail(pca_data)
dim(pca_data)
# Add group information to PCA data
pca_data$group <- c(rep("WT", 9), rep("KO_AB16", 9),
                    rep("KO_AE4", 9), rep("NC_ctrl", 8))

# Plot PCA
# pretty label - sub("^[^_]+_", "", rownames
ggplot(pca_data, aes(x = PC1, y = PC2, color = group, label = rownames(pca_data))) +
  geom_point() +
  labs(title = "PCA, 3 technical reps * 3 biological reps for each group")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_text_repel(size = 3)

# clustering
rnames <- c("Parental_1.1", "Parental_1.2", "Parental_1.3",
            "Parental_2.1", "Parental_2.2", "Parental_2.3",
            "Parental_3.1", "Parental_3.2", "Parental_3.3",
            "AB16_1.1", "AB16_1.2", "AB16_1.3", "AB16_2.1", "AB16_2.2", "AB16_2.3",
            "AB16_3.1", "AB16_3.2", "AB16_3.3",
            "AE4_1.1", "AE4_1.2", "AE4_1.3", "AE4_2.1", "AE4_2.2", "AE4_2.3",
            "AE4_3.1", "AE4_3.2", "AE4_3.3",
            "NC2_1.1", "NC2_1.2", "NC2_1.3", "NC2_2.1", "NC2_2.2", "NC2_2.3",
            "NC2_3.1", "NC2_3.2", "NC2_3.3")

samples <- data.frame(sample = as.character(colnames(intens[,12:38])), Sample_type="METTL1_KO") %>%
  column_to_rownames("sample")
samples$Sample_type <- c(#rep("Parental_1", 3), rep("Parental_2", 3), rep("Parental_3", 3), 
                         rep("AB16_1", 3), rep("AB16_2", 3), rep("AB16_3", 3),
                         rep("AE4_1", 3), rep("AE4_2", 3), rep("AE4_3", 3), 
                         rep("NC2_1", 3), rep("NC2_2", 3), rep("NC2_3", 3))
cols <- list(Sample_type=c(#Parental_1="#FFD1D1",Parental_2="#FFB3B3", Parental_3="#FF9696",
                AB16_1="#CCE5FF", AB16_2 ="#99C2FF", AB16_3="#66A0FF",
                AE4_1="#C9FFD1", AE4_2 ="#9DFFA6", AE4_3="#71FF7B", 
                NC2_1="#D2B48C", NC2_2="#C19A6B", NC2_3="#A97C4D"))

sampleDists <- dist(t(intens[,12:38]))
sampleDistMatrix <- as.matrix( sampleDists )
heatmap <- pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists,
               show_rownames = F, show_colnames = F,
               annotation_row = samples, annotation_col = samples,
         annotation_colors = cols,annotation_names_col=F,annotation_names_row=F,legend=F,
         main = "Proteomics samples similarities")

############### selecting background difference - nc_wt vs parental 

ctrlp_ctrl_nc <- filter(ratio, (compare_to == "OE33_Parental_WT") & (compare_from == "NC2_WT"))
# 19 in total:
# ctrlp_ctrl_nc_signif <- filter(ctrlp_ctrl_nc, (adj_pvalue < 0.05) & (abs(log2FC) > 2)) 
# 159 in total:
ctrlp_ctrl_nc_signif <- filter(ctrlp_ctrl_nc, (adj_pvalue < 0.05) & (abs(log2FC) > 1))
# 865 in total:
#ctrlp_ctrl_nc_signif <- filter(ctrlp_ctrl_nc, (adj_pvalue < 0.05)) 
head(ctrlp_ctrl_nc_signif)
dim(ctrlp_ctrl_nc_signif)

# creating heatmap for ctrl-ctrl comparison
ctrlp_ctrl_nc_signif <- merge(x = ctrlp_ctrl_nc_signif, 
                     y = intens[,c(2:11, 30:38)], 
                     by.x="Symbol",
                     by.y="symbol")

# removing the proteins where all the values = 0 
ctrlp_ctrl_nc_signif <- ctrlp_ctrl_nc_signif %>%
  filter(rowSums(.[, 11:28]) != 0)

samples <- data.frame(sample = as.character(colnames(ctrlp_ctrl_nc_signif[, 11:28])),
                      Sample_type="Exp") %>% column_to_rownames("sample")
samples$Sample_type <- c(rep("parental_wt", 9), rep("nc_wt", 9))

hm_gg <- as.ggplot(pheatmap(as.matrix(log(ctrlp_ctrl_nc_signif[, 11:28]+1)), scale="row",
                            cluster_rows = TRUE,
                            fontsize_col = 8,
                            angle_col = 90,
                            cellheight = 2,
                            treeheight_col = 3,
                            treeheight_row = 0,
                            annotation_col=samples,
                            annotation_colors =list(
                              Sample_type=c(parental_wt="orange",nc_wt="black")),
                            color=colorRampPalette(c("navy","white", "red"))(50), 
                            show_rownames = F, annotation_names_col = FALSE))

# proteins list TO SUBSTRACT FROM OTHER COMPS
prots_background <- ctrlp_ctrl_nc_signif$Symbol
length(prots_background)

######################## selecting each comparison 

# parsing raw stats for comparisons

unique(raw_stats$Label)
# KO to parental
comps <- raw_stats %>% select("Protein", "Label", "log2FC", "adj.pvalue") %>% filter(
  Label %in% c("AB16_METTL1_KO-OE33_Parental_WT", "AE4_METTL1_KO-OE33_Parental_WT", 
               "AE4_METTL1_KO-NC2_WT", "AB16_METTL1_KO-NC2_WT")
  )
comps <- comps %>%
  separate(Label, into = c("compare_to", "compare_from"), sep = "-")
comps$Symbol <- sub("^.*\\|([^_]+)_HUMAN$", "\\1", comps$Protein)
head(comps$Protein)
ab16_ctrlp <- filter(comps, (compare_to == "AB16_METTL1_KO") & (compare_from == "OE33_Parental_WT"))
ae4_ctrlp <- filter(comps, (compare_to == "AE4_METTL1_KO") & (compare_from == "OE33_Parental_WT"))
ab16_ctrl_nc <- filter(comps, (compare_to == "AB16_METTL1_KO") & (compare_from == "NC2_WT"))
ae4_ctrl_nc <- filter(comps, (compare_to == "AE4_METTL1_KO") & (compare_from == "NC2_WT"))
head(ab16_ctrlp)

# KO to nc - from ratio
#ab16_ctrl_nc <- filter(ratio, (compare_to == "AB16_METTL1_KO") & (compare_from == "NC2_WT"))
#ae4_ctrl_nc <- filter(ratio, (compare_to == "AE4_METTL1_KO") & (compare_from == "NC2_WT"))
#head(ab16_ctrl_nc)


########### Functions for heatmap ###########

create_heatmap <- function(dataset, col_numbers_in_iner, gr_names){
  dt_signif <- filter(dataset,
                      (adj.pvalue < 0.05) & (abs(log2FC) != Inf))
  #head(dt_signif)
  #dim(dt_signif)
  
  # creating heatmap for ctrl-ctrl comparison
  dt_signif <- merge(x = dt_signif, y = intens[,col_numbers_in_iner], 
                     by.x="Symbol", by.y="symbol")
  
  # removing the proteins where all the values = 0 and fc != 0
  dt_signif <- dt_signif %>% filter((rowSums(.[, 7:24]) != 0))
  dt_signif_sorted <- arrange(dt_signif, desc(log2FC))
  top_200_rows <- head(dt_signif_sorted, n = min(nrow(dt_signif_sorted), 150))
  
  samples <- data.frame(sample = as.character(colnames(top_200_rows[, 7:24])),
                        Sample_type="Exp") %>% column_to_rownames("sample")
  samples$Sample_type <- c(rep(gr_names[1], 9), rep(gr_names[2], 9))
  
  hm_gg <- as.ggplot(pheatmap(as.matrix(log(top_200_rows[, 7:24]+1)), scale="row",
                              cluster_rows = TRUE,
                              fontsize_col = 8,
                              angle_col = 90,
                              cellheight = 2,
                              treeheight_col = 3,
                              treeheight_row = 0,
                              annotation_col=samples,
                              #annotation_colors =list(
                              #  Sample_type=c(gr_names[1]="orange",gr_names[2]="black")),
                              color=colorRampPalette(c("navy","white", "red"))(50), 
                              show_rownames = F, annotation_names_col = FALSE))
  return(hm_gg)
}

############  AB16 - parent ctrl #####

ab16_ctrlp_for_inters <- filter(ab16_ctrlp, (adj.pvalue < 0.05))
ab16_ctrlp_for_inters_pr <- setdiff(ab16_ctrlp_for_inters$Symbol, prots_background)
dim(ab16_ctrlp_for_inters)
length(ab16_ctrlp_for_inters_pr)
# from 892 to 867

hm <- create_heatmap(ab16_ctrlp, c(2:20), c("AB16", "parental_wt"))

############  AE4 - parent ctrl #####

ae4_ctrlp_for_inters <- filter(ae4_ctrlp, (adj.pvalue < 0.05))
ae4_ctrlp_for_inters_pr <- setdiff(ae4_ctrlp_for_inters$Symbol, prots_background)
dim(ae4_ctrlp_for_inters)
length(ae4_ctrlp_for_inters_pr)
# from 1201 to 1165

hm <- create_heatmap(ae4_ctrlp, c(2:11, 21:30), c("AE4", "parental_wt"))

############  AB16 - nc wt #####

ab16_nc_for_inters <- filter(ab16_ctrl_nc, (adj.pvalue < 0.05))
ab16_nc_for_inters_pr <- setdiff(ab16_nc_for_inters$Symbol, prots_background)
dim(ab16_nc_for_inters)
length(ab16_nc_for_inters_pr)
# from 1247 to 1210

hm <- create_heatmap(ab16_ctrl_nc, c(2, 12:20, 29:38), c("AB16", "NC_WT"))

############  AE4 - nc wt #####

ae4_nc_for_inters <- filter(ae4_ctrl_nc, (adj.pvalue < 0.05))
ae4_nc_for_inters_pr <- setdiff(ae4_nc_for_inters$Symbol, prots_background)
dim(ae4_nc_for_inters)
length(ae4_nc_for_inters_pr)
# from 1046 to 1020

hm <- create_heatmap(ae4_ctrl_nc, c(2, 20:38), c("AE4", "NC_WT"))









############ Overlaps with key proteins WO parental ##########
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

############ Overlaps with key proteins ONLY 3RD REPS ##########

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


















############ CORRECT DE analysis for for proteomics data ##############

# using only 3rd reps and w/o parental 
# we will do two approaches:
# 1) select intersection bw AB16-NC2 & AE4-NC2; build network --> key prots
# 2) build 2 networks for each comparison (AB16-NC2, AE4-NC2) --> get 2 key prots sets --> intersect them

# DATA PREP ----------------------------------------------------------------

# selecting cols 
intens_for_kp <-  intens[,c(2, 18:20, 27:29, 36:38)]
head(intens_for_kp)
dim(intens_for_kp) # 3932 10
# removing rows with all 0s
intens_for_kp <- filter(intens_for_kp, rowSums(intens_for_kp[-1]) != 0)
dim(intens_for_kp) # 3633 10
# renaming cols 
colnames(intens_for_kp) <- c("symbol","AB16_METTL1_KO_3.1","AB16_METTL1_KO_3.2","AB16_METTL1_KO_3.3",
                             "AE4_METTL1_KO_3.1","AE4_METTL1_KO_3.2","AE4_METTL1_KO_3.3",
                             "NC2_WT_3.1","NC2_WT_3.2","NC2_WT_3.3")

# using SummirizedExperiment to store vals
# based on https://rformassspectrometry.github.io/book/sec-quant.html 

symbols <- intens_for_kp[, 1] # First column contains peptide names
pep_counts <- as.matrix(intens_for_kp[, -1]) # Exclude the first column for counts
# creating QFeatures data
cptac <- readQFeatures(intens_for_kp, ecol = 2:10, name = "proteins", fnames = "symbol")
cptac$group <- c(rep("AB16_METTL1_KO", 3),
                  rep("AE4_METTL1_KO", 3), rep("NC2_WT", 3))
colData(cptac)
rowDataNames(cptac)
head(assay(cptac))

# zeros as NA -> WE ALREADY REMOVED ALL = 0
cptac <- zeroIsNA(cptac, i = seq_along(cptac))
nNA(cptac, i = seq_along(cptac))
#barplot(nNA(cptac, i = seq_along(cptac))$nNAcols$pNA)
table(nNA(cptac, i = seq_along(cptac))$nNArows$nNA)
# remove rows that have >=7 NAs out of 9
cptac <- filterNA(cptac, i = seq_along(cptac), pNA = 7/9) # 3374 prots left

# trying IMPUTATION ----------------------------------------------------------

# we have 700+ prots with NAs
cls <- c("black", "red", "blue", "steelblue", "orange")
plot(density(na.omit(assay(cptac))), col = cls[1],  xlim = c(0, 4e+09), ylim = c(0, 5.0e-09))
lines(density(assay(impute(cptac, method = "knn"))), col = cls[2])
lines(density(assay(impute(cptac, method = "zero"))), col = cls[3])
lines(density(assay(impute(cptac, method = "MinDet"))), col = cls[4])
lines(density(assay(impute(cptac, method = "bpca"))), col = cls[5])
legend("topright", legend = c("orig", "knn", "zero", "MinDet", "bpca"),
       col = cls, lwd = 2, bty = "n")
# from worst to best: bpca, MinDet = zero, knn, orig 
cls <- c("black", "red")
plot(density(na.omit(assay(cptac))), col = cls[1],  xlim = c(0, 2e+09), ylim = c(0, 4.0e-09))
lines(density(assay(impute(cptac, method = "knn"))), col = cls[2])
legend("topright", legend = c("orig", "knn"),
       col = cls, lwd = 2, bty = "n")

# WE WILL PERFORM TWO SIMULTANEOUS ANALYSIS - WITH MinDet AND W/O IMPUTATION
# MinDet looks the best but gives terrible distribution after log norm with a
# 2nd peak in lower values; same for min imputation.
# knn and bpca don't shift the distribution!
cptac_imp <- impute(cptac, method = "bpca")
cptac <- impute(cptac, method = "zero")

# counting unique features
cptac <- countUniqueFeatures(cptac,
                             i = "proteins",
                             colDataName = "proteins_counts")
colData(cptac)
??countUniqueFeatures

# TRANSFORMATIONS -----------------------------------------------------------

# log-transformation
cptac <- addAssay(cptac,
                  logTransform(cptac[[1]]),
                  name = "prots_log")
cptac_imp <- addAssay(cptac_imp,
                  logTransform(cptac_imp[[1]]),
                  name = "prots_log")

# normalization
cptac <- addAssay(cptac,
                  normalize(cptac[["prots_log"]], method = "center.median"),
                  name = "prots_norm")
cptac_imp <- addAssay(cptac_imp,
                  normalize(cptac_imp[["prots_log"]], method = "center.median"),
                  name = "prots_norm")

# Plotting
par(mfrow = c(1, 3))
limma::plotDensities(assay(cptac[["proteins"]]))
limma::plotDensities(assay(cptac[["prots_log"]]))
limma::plotDensities(assay(cptac[["prots_norm"]]))

par(mfrow = c(1, 3))
limma::plotDensities(assay(cptac_imp[["proteins"]]))
limma::plotDensities(assay(cptac_imp[["prots_log"]]))
limma::plotDensities(assay(cptac_imp[["prots_norm"]]))

# DE ITSELF ------------------------------------------------------------------

# using limma to calculate p-vals and lfc
# limma will output NA in p-cal if protein has NA is more that 2 groups 

prots <- getWithColData(cptac, "prots_norm")
prots_imputed <- getWithColData(cptac_imp, "prots_norm")

groups <- prots$group
design <- model.matrix(~ 0 + groups)
contrast_matrix <- makeContrasts(
  AB16_METTL1_KO_vs_NC2_WT = groupsAB16_METTL1_KO - groupsNC2_WT,
  AE4_METTL1_KO_vs_NC2_WT = groupsAE4_METTL1_KO - groupsNC2_WT,
  levels = colnames(design)
)

# zeros
fit_prot <- lmFit(assay(prots), design)
contrast_fit_prot <- contrasts.fit(fit_prot, contrast_matrix)
eb_fit_prot <- eBayes(contrast_fit_prot)
# imputed
fit_prot_imp <- lmFit(assay(prots_imputed), design)
contrast_fit_prot_imp <- contrasts.fit(fit_prot_imp, contrast_matrix)
eb_fit_prot_imp <- eBayes(contrast_fit_prot_imp)

# zeros
res_AB16 <-
  topTable(eb_fit_prot, coef = "AB16_METTL1_KO_vs_NC2_WT", number = Inf) %>%
  rownames_to_column("protein") %>%
  as_tibble() %>%
  mutate(TP = grepl("ups", protein))
res_AE4 <-
  topTable(eb_fit_prot, coef = "AE4_METTL1_KO_vs_NC2_WT", number = Inf) %>%
  rownames_to_column("protein") %>%
  as_tibble() %>%
  mutate(TP = grepl("ups", protein))
# imputed
res_AB16_imp <-
  topTable(eb_fit_prot_imp, coef = "AB16_METTL1_KO_vs_NC2_WT", number = Inf) %>%
  rownames_to_column("protein") %>%
  as_tibble() %>%
  mutate(TP = grepl("ups", protein))
res_AE4_imp <-
  topTable(eb_fit_prot_imp, coef = "AE4_METTL1_KO_vs_NC2_WT", number = Inf) %>%
  rownames_to_column("protein") %>%
  as_tibble() %>%
  mutate(TP = grepl("ups", protein))

# SAVING THE RESULTS ----------------------------------------------------------

# saving the full data
write.csv(res_AE4, file = "./reps3nc2_analysis/full_data/res_AE4.csv", row.names = FALSE)
write.csv(res_AE4_imp, file = "./reps3nc2_analysis/full_data/res_AE4_imp.csv", row.names = FALSE)
write.csv(res_AB16_imp, file = "./reps3nc2_analysis/full_data/res_AB16_imp.csv", row.names = FALSE)
write.csv(res_AB16, file = "./reps3nc2_analysis/full_data/res_AB16.csv", row.names = FALSE)

# filtering the data
res_AE4_imp <- res_AE4_imp %>% filter(adj.P.Val < 0.05 & abs(logFC)>0.5)
res_AE4 <- res_AE4 %>% filter(adj.P.Val < 0.05 & abs(logFC)>0.5)
res_AB16_imp <- res_AB16_imp %>% filter(adj.P.Val < 0.05 & abs(logFC)>0.5)
res_AB16 <- res_AB16 %>% filter(adj.P.Val < 0.05 & abs(logFC)>0.5)

# saving the filtered data
write.csv(res_AE4, file = "./reps3nc2_analysis/filtered/res_AE4.csv", row.names = FALSE)
write.csv(res_AE4_imp, file = "./reps3nc2_analysis/filtered/res_AE4_imp.csv", row.names = FALSE)
write.csv(res_AB16_imp, file = "./reps3nc2_analysis/filtered/res_AB16_imp.csv", row.names = FALSE)
write.csv(res_AB16, file = "./reps3nc2_analysis/filtered/res_AB16.csv", row.names = FALSE)
write.table(res_AE4$protein, file = "./reps3nc2_analysis/filtered/res_AE4_names.txt",
            row.names = FALSE, quote = F, col.names = F)
write.table(res_AE4_imp$protein, file = "./reps3nc2_analysis/filtered/res_AE4_imp_names.txt",
            row.names = FALSE, quote = F, col.names = F)
write.table(res_AB16$protein, file = "./reps3nc2_analysis/filtered/res_AB16_names.txt",
            row.names = FALSE, quote = F, col.names = F)
write.table(res_AB16_imp$protein, file = "./reps3nc2_analysis/filtered/res_AB16_imp_names.txt",
            row.names = FALSE, quote = F, col.names = F)

# INTERSECTIONS between comparisons -------------------------------------------

# samples intersection
inters_AB16_AE4_zeros_names <- intersect(res_AB16$protein, res_AE4$protein)
inters_AB16_AE4_imp_names <- intersect(res_AB16_imp$protein, res_AE4_imp$protein)

# adding col with mean logFC
calc_mean_lfc_prots <- function(prnames, source1, source2){
  new_data <- data.frame(prnames)
  new_data$mean_lfc <- sapply(new_data$prnames, function(name) {
    mean(c(subset(source1, source1$protein==name)$logFC,
           subset(source2, source2$protein==name)$logFC))
  })
  new_data$same_direction <- sapply(new_data$prnames, function(name){
    sign(subset(source1, source1$protein==name)$logFC) == sign(subset(source2, source2$protein==name)$logFC)
  })
  return(new_data)
}

inters_AB16_AE4_zeros_vals <- calc_mean_lfc_prots(inters_AB16_AE4_zeros_names, res_AB16, res_AE4)
inters_AB16_AE4_imp_vals <- calc_mean_lfc_prots(inters_AB16_AE4_imp_names, res_AB16_imp, res_AE4_imp)

# saving
write.csv(inters_AB16_AE4_zeros_vals,
          file = "./reps3nc2_analysis/filtered/intersections/inters_AB16_AE4_zeros_vals.csv",
          row.names = FALSE)
write.csv(inters_AB16_AE4_imp_vals,
          file = "./reps3nc2_analysis/filtered/intersections/inters_AB16_AE4_imp_vals.csv",
          row.names = FALSE)
write.table(inters_AB16_AE4_zeros_names,
            file = "./reps3nc2_analysis/filtered/intersections/inters_AB16_AE4_zeros_names.txt",
            row.names = FALSE, quote = F, col.names = F)
write.table(inters_AB16_AE4_imp_names,
            file = "./reps3nc2_analysis/filtered/intersections/inters_AB16_AE4_imp_names.txt",
            row.names = FALSE, quote = F, col.names = F)

# plotting the venn diagrams for intersections
gene_sets <- list(AB16_NC2 = res_AB16$protein, AE4_NC2 = res_AE4$protein)
circle_colors <- c("AB16_NC2" = "skyblue", "AE4_NC2" = "lightcoral")
euler_result <- euler(gene_sets, shape = "ellipse", fills = circle_colors)
plot(euler_result, quantities = TRUE, fills = circle_colors, counts = TRUE,
     legend = TRUE, main = "Not imputed: AB16-NC2 vs AE4-NC2")

gene_sets <- list(AB16_NC2 = res_AB16_imp$protein, AE4_NC2 =  res_AE4_imp$protein)
circle_colors <- c("AB16_NC2" = "skyblue", "AE4_NC2" = "lightcoral")
euler_result <- euler(gene_sets, shape = "ellipse", fills = circle_colors)
plot(euler_result, quantities = TRUE, fills = circle_colors, counts = TRUE,
     legend = TRUE, main = "Imputed: AB16-NC2 vs AE4-NC2")

gene_sets <- list(AE4_NC2_imp = res_AE4_imp$protein, AE4_NC2 =  res_AE4$protein)
circle_colors <- c("AE4_NC2_imp" = "skyblue", "AE4_NC2" = "lightcoral")
euler_result <- euler(gene_sets, shape = "ellipse", fills = circle_colors)
plot(euler_result, quantities = TRUE, fills = circle_colors, counts = TRUE,
     legend = TRUE, main = "Imputed vs Not-Imp, AE4")

gene_sets <- list(AB16_NC2_imp = res_AB16_imp$protein, AB16_NC2 =  res_AB16$protein)
circle_colors <- c("AB16_NC2_imp" = "skyblue", "AB16_NC2" = "lightcoral")
euler_result <- euler(gene_sets, shape = "ellipse", fills = circle_colors)
plot(euler_result, quantities = TRUE, fills = circle_colors, counts = TRUE,
     legend = TRUE, main = "Imputed vs Not-Imp, AB16")

################# Full RNA-seq - full prots intersection ################

# raw data inters
full_rnaseq_prots_names_t <- intersect(contrast1t$genes$gene_name, intens_for_kp$symbol) # 3329
full_rnaseq_prots_names_g <- intersect(contrast1g$genes$gene_name, intens_for_kp$symbol) # 3331
length(full_rnaseq_prots_names_g)

# intersections with the comparisons filtered by p-val<0.05 and lfc>0.5, INT OF EDGER and DESEQ2
ae4_ab1_gpinters <- intersect(ae4_ab1_intg_full$gene_name, intens_for_kp$symbol) # 886
ab16_ab1_gpinters <- intersect(ab16_ab1_intg_full$gene_name, intens_for_kp$symbol) # 571
ab16_ae21_gpinters <-intersect(ab16_ae21_intg_full$gene_name, intens_for_kp$symbol) # 760
ae4_ae21_gpinters <-intersect(ae4_ae21_intg_full$gene_name, intens_for_kp$symbol) # 1179

ae4_ab1_tpinters <-intersect(ae4_ab1_intt_full$gene_name, intens_for_kp$symbol) # 427
ab16_ab1_tpinters <-intersect(ab16_ab1_intt_full$gene_name, intens_for_kp$symbol) # 279
ab16_ae21_tpinters <-intersect(ab16_ae21_intt_full$gene_name, intens_for_kp$symbol) # 349
ae4_ae21_tpinters <-intersect(ae4_ae21_intt_full$gene_name, intens_for_kp$symbol) # 504

gi_g_pinters <-intersect(global_int_g_counts$gene_name, intens_for_kp$symbol) # 114
gi_t_pinters <-intersect(global_int_t_counts$gene_name, intens_for_kp$symbol) # 34

# union
genes_union <- union(union(ae4_ab1_gpinters, ab16_ab1_gpinters), union(ab16_ae21_gpinters, ae4_ae21_gpinters))
length(genes_union) # 1821
trs_union <- union(union(ae4_ab1_tpinters, ab16_ab1_tpinters), union(ab16_ae21_tpinters, ae4_ae21_tpinters))
length(trs_union) # 946

# same direction rna-seq
down_genes <- select_one_direction(genes_union)[['down']] # 598
up_genes <- select_one_direction(genes_union)[['up']] # 1167
down_trs <- select_one_direction(trs_union)[['down']] # 283
up_trs <- select_one_direction(trs_union)[['up']] # 633
genes_union <- append(down_genes, up_genes)
trs_union <- append(down_trs, up_trs)
# same direction proteomics
dfg <- data.frame(genes_union)
dft <- data.frame(trs_union)

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
dim(filter(dfg, same_direction==TRUE))
dim(filter(dfg, same_direction_imp==TRUE))
dim(filter(dft, same_direction==TRUE))
dim(filter(dft, same_direction_imp==TRUE))

# saving to df the last step of the early intersection
df_fin_t <- filter(dft, same_direction_imp==TRUE)
dim(df_fin_t) # 557
df_fin_g <- filter(dfg, same_direction_imp==TRUE)
dim(df_fin_g) # 1055

# saving for the networks construction
write.table(df_fin_t$trs_union,
            file = "early_prot_inters_trs_names.txt",
            row.names = FALSE, quote = F, col.names = F)
write.table(df_fin_g$genes_union,
            file = "early_prot_inters_genes_names.txt",
            row.names = FALSE, quote = F, col.names = F)

# adding mean lcf for spia
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

# other intersections

ky_prots_proteomics <- read.table("./reps3nc2_analysis/key_prots_total.txt", header = FALSE)$V1
int<- intersect(global_int_t_genenames, inters_AB16_AE4_imp_names)
init_tr_p_inters <- c('CAD', 'DTYMK', 'FEN1', 'MCM3', 'MCM4', 'MCM6', 'MAT2A',
                      'MRPL4', 'NOP56', 'NCL', 'POLR1A', 'ZWINT')
print(filter(inters_AB16_AE4_imp_vals, prnames %in% init_tr_p_inters))
head(inters_AB16_AE4_imp_vals)

# intersection between transcr, prot, and genes 

three_inters <- data.frame()
names_in_both <- intersect(df_fin_t$trs_union, df_fin_g$genes_union)
length(names_in_both) # 479

for(i in 1:nrow(df_fin_t)) {
  if((df_fin_t$trs_union[i] %in% names_in_both) && (sign(df_fin_t$mean_lfc[i]) == sign(df_fin_g$mean_lfc[i]))){
    row_to_add <- data.frame(gene_name = df_fin_t$trs_union[i], mean_lfc_tr = df_fin_t$mean_lfc[i], mean_lfc_g = df_fin_g$mean_lfc[i])
    three_inters <- rbind(three_inters, row_to_add)
  }
}

head(three_inters)
dim(three_inters) # 240

write.table(three_inters$gene_name,
            file = "early_prot_3inters_names.txt",
            row.names = FALSE, quote = F, col.names = F)







