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

########################## Preliminary data exploration ########################

# In order to explore the distribution of the intensity counts and any differences
# between samples, we will plot these distributions on the log-scale: 
# log-scaled intensity against number of peptides with this intensity value.

# Saving non-zero minimal values for each column of the intensity data (intens)
# to the non0mins vector. We will later use this vector to plot these values on
# the distribution plot.
non0mins <- c()
for(col in colnames(intens[,3:38])){
  non0mins <- append(non0mins, log(unique(sort(intens[,col]))[2]))
}
non0mins <- as.data.frame(non0mins)
rownames(non0mins) <- colnames(intens[,3:38])
colnames(non0mins) <- c('value')

# Normalizing the intensity counts - log transform
long_f_log_wt <- gather(log(as.data.frame(intens[,3:11] +1)))
long_f_log_ko_ab16 <- gather(log(as.data.frame(intens[,12:20] +1)))
long_f_log_ko_ae4 <- gather(log(as.data.frame(intens[,21:29] +1)))
long_f_log_nc2_wt <- gather(log(as.data.frame(intens[,30:38] +1)))
head(intens[,30:38])
dim(long_f_log_wt)

# Function to add min values to the intens dataframes
create_vec <- function(vec, n){
  cout <- c()
  for(val in vec){
    cout <- append(cout, rep(val, n))
  }
  return(cout)
}

# Adding non-zero min values to intensity dataframes 
long_f_log_wt$min <- create_vec(non0mins[1:9,], 3932)
long_f_log_ko_ab16$min <- create_vec(non0mins[10:18,], 3932)
long_f_log_ko_ae4$min <- create_vec(non0mins[19:27,], 3932)
long_f_log_nc2_wt$min <- create_vec(non0mins[28:36,], 3932)

# Plotting the log-normalized intensities for different samples
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


################################# PCA ##########################################

# Deleting rows with all zeros, +1 because the sample X36 is an outlier
all_log <- log(intens[,3:37]+1)
all_log <- subset(all_log, rowSums(all_log != 0) > 0)

# Computing the PCA
pc <- prcomp(t(as.matrix(all_log)), center = TRUE, scale = TRUE)
pca_data <- as.data.frame(pc$x)
# Adding group information to PCA data
pca_data$group <- c(rep("WT", 9), rep("KO_AB16", 9),
                    rep("KO_AE4", 9), rep("NC_ctrl", 8))

# Plot PCA
# to make the labels prettier: sub("^[^_]+_", "", rownames)
ggplot(pca_data, aes(x = PC1, y = PC2, color = group, label = rownames(pca_data))) +
  geom_point() +
  labs(title = "PCA, 3 technical reps * 3 biological reps for each group")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_text_repel(size = 3)

############################## Samples clustering ##############################

# Below we are computing the distance matrix between samples, and plotting it 
# along with hierarchical clustering of the samples.

# Renaming samples for clustering
rnames <- c("Parental_1.1", "Parental_1.2", "Parental_1.3",
            "Parental_2.1", "Parental_2.2", "Parental_2.3",
            "Parental_3.1", "Parental_3.2", "Parental_3.3",
            "AB16_1.1", "AB16_1.2", "AB16_1.3", "AB16_2.1", "AB16_2.2", "AB16_2.3",
            "AB16_3.1", "AB16_3.2", "AB16_3.3",
            "AE4_1.1", "AE4_1.2", "AE4_1.3", "AE4_2.1", "AE4_2.2", "AE4_2.3",
            "AE4_3.1", "AE4_3.2", "AE4_3.3",
            "NC2_1.1", "NC2_1.2", "NC2_1.3", "NC2_2.1", "NC2_2.2", "NC2_2.3",
            "NC2_3.1", "NC2_3.2", "NC2_3.3")

# Data preparation
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

# Computing the distance matrix
sampleDists <- dist(t(intens[,12:38]))
sampleDistMatrix <- as.matrix( sampleDists )

# Plotting
heatmap <- pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists,
               show_rownames = F, show_colnames = F,
               annotation_row = samples, annotation_col = samples,
         annotation_colors = cols,annotation_names_col=F,annotation_names_row=F,legend=F,
         main = "Proteomics samples similarities")

###################### Heatmaps for protein intensities ######################## 

####### NC2 - Parental heatmap

# First, we will look closer at the differences between two controls:
# NC2_wt versus parental samples. Generally, the difference between these
# controls should serve as a background noise/batch difference that we would
# need to subtract from all other METTL1 KO vs controls comparisons. 

ctrlp_ctrl_nc <- filter(ratio, (compare_to == "OE33_Parental_WT") & (compare_from == "NC2_WT"))

# In order to select proteins to plot we can apply p-adjusted and logFC cut-offs:

# * 19 proteins in total with 0.05 and 2 respectfully
# ctrlp_ctrl_nc_signif <- filter(ctrlp_ctrl_nc, (adj_pvalue < 0.05) & (abs(log2FC) > 2)) 

# * 159 proteins in total with 0.05 and 1 respectfully
ctrlp_ctrl_nc_signif <- filter(ctrlp_ctrl_nc, (adj_pvalue < 0.05) & (abs(log2FC) > 1))

# * 865 proteins in total with only p-adjusted < 0.05
#ctrlp_ctrl_nc_signif <- filter(ctrlp_ctrl_nc, (adj_pvalue < 0.05))


# Creating a heatmap for this NC2-WT vs Parental comparison
ctrlp_ctrl_nc_signif <- merge(x = ctrlp_ctrl_nc_signif, 
                     y = intens[,c(2:11, 30:38)], 
                     by.x="Symbol",
                     by.y="symbol")

# Removing the proteins where all the values = 0 
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


####### Saving each comparison to separated files 

# Parsing raw stats for comparisons
# KO to parental
comps <- raw_stats %>% select("Protein", "Label", "log2FC", "adj.pvalue") %>% filter(
  Label %in% c("AB16_METTL1_KO-OE33_Parental_WT", "AE4_METTL1_KO-OE33_Parental_WT", 
               "AE4_METTL1_KO-NC2_WT", "AB16_METTL1_KO-NC2_WT")
  )
comps <- comps %>%
  separate(Label, into = c("compare_to", "compare_from"), sep = "-")
comps$Symbol <- sub("^.*\\|([^_]+)_HUMAN$", "\\1", comps$Protein)

ab16_ctrlp <- filter(comps, (compare_to == "AB16_METTL1_KO") & (compare_from == "OE33_Parental_WT"))
ae4_ctrlp <- filter(comps, (compare_to == "AE4_METTL1_KO") & (compare_from == "OE33_Parental_WT"))
ab16_ctrl_nc <- filter(comps, (compare_to == "AB16_METTL1_KO") & (compare_from == "NC2_WT"))
ae4_ctrl_nc <- filter(comps, (compare_to == "AE4_METTL1_KO") & (compare_from == "NC2_WT"))

# KO to NC2 - from ratio dataset
#ab16_ctrl_nc <- filter(ratio, (compare_to == "AB16_METTL1_KO") & (compare_from == "NC2_WT"))
#ae4_ctrl_nc <- filter(ratio, (compare_to == "AE4_METTL1_KO") & (compare_from == "NC2_WT"))


####### Function for the heatmap

# This function takes the dataset with intensities, column numbers representing 
# the samples we want to compare, and the names of the groups we are comparing,
# and creates a heatmap plot. 

create_heatmap <- function(dataset, col_numbers_in_iner, gr_names){
  dt_signif <- filter(dataset,
                      (adj.pvalue < 0.05) & (abs(log2FC) != Inf))
  
  # Creating heatmap for ctrl-ctrl comparison
  dt_signif <- merge(x = dt_signif, y = intens[,col_numbers_in_iner], 
                     by.x="Symbol", by.y="symbol")
  
  # Removing the proteins where all the values = 0 and fc != 0
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

####### AB16 - parent ctrl

ab16_ctrlp_for_inters <- filter(ab16_ctrlp, (adj.pvalue < 0.05))
ab16_ctrlp_for_inters_pr <- setdiff(ab16_ctrlp_for_inters$Symbol, prots_background)
#dim(ab16_ctrlp_for_inters)
#length(ab16_ctrlp_for_inters_pr)
# from 892 to 867

hm <- create_heatmap(ab16_ctrlp, c(2:20), c("AB16", "parental_wt"))

####### AE4 - parent ctrl

ae4_ctrlp_for_inters <- filter(ae4_ctrlp, (adj.pvalue < 0.05))
ae4_ctrlp_for_inters_pr <- setdiff(ae4_ctrlp_for_inters$Symbol, prots_background)
#dim(ae4_ctrlp_for_inters)
#length(ae4_ctrlp_for_inters_pr)
# from 1201 to 1165

hm <- create_heatmap(ae4_ctrlp, c(2:11, 21:30), c("AE4", "parental_wt"))

####### AB16 - NC2 wt

ab16_nc_for_inters <- filter(ab16_ctrl_nc, (adj.pvalue < 0.05))
ab16_nc_for_inters_pr <- setdiff(ab16_nc_for_inters$Symbol, prots_background)
#dim(ab16_nc_for_inters)
#length(ab16_nc_for_inters_pr)
# from 1247 to 1210

hm <- create_heatmap(ab16_ctrl_nc, c(2, 12:20, 29:38), c("AB16", "NC_WT"))

####### AE4 - NC2 wt

ae4_nc_for_inters <- filter(ae4_ctrl_nc, (adj.pvalue < 0.05))
ae4_nc_for_inters_pr <- setdiff(ae4_nc_for_inters$Symbol, prots_background)
#dim(ae4_nc_for_inters)
#length(ae4_nc_for_inters_pr)
# from 1046 to 1020

hm <- create_heatmap(ae4_ctrl_nc, c(2, 20:38), c("AE4", "NC_WT"))



########## Differential expression analysis for for proteomics data ############

# As we saw on the correlation map and the PCA plot, samples are quire heterogeneous 
# this to catch the reasonable protein level changes between knockout samples and 
# controls, we will be using only 3rd repetisions and only NC2_wt control.

# We will perofrm the analysis in two ways:
# 1) select intersection between AB16-NC2 & AE4-NC2
#    --> build network --> extract key proteins
# 2) build two networks for each comparison (AB16-NC2, AE4-NC2)
#    --> extract two key protein sets --> intersect them


# DATA PREPARATION -------------------------------------------------------------

# Selecting columns 
intens_for_kp <-  intens[,c(2, 18:20, 27:29, 36:38)]
head(intens_for_kp)
dim(intens_for_kp) # 3932 10
# Removing rows with all zeros
intens_for_kp <- filter(intens_for_kp, rowSums(intens_for_kp[-1]) != 0)
dim(intens_for_kp) # 3633 10
# Renaming columns
colnames(intens_for_kp) <- c("symbol","AB16_METTL1_KO_3.1","AB16_METTL1_KO_3.2","AB16_METTL1_KO_3.3",
                             "AE4_METTL1_KO_3.1","AE4_METTL1_KO_3.2","AE4_METTL1_KO_3.3",
                             "NC2_WT_3.1","NC2_WT_3.2","NC2_WT_3.3")

# Here we are using sing QFeatures objects to store values
# based on https://rformassspectrometry.github.io/book/sec-quant.html 

# First column contains peptide names
symbols <- intens_for_kp[, 1] 
# Exclude the first column for counts
pep_counts <- as.matrix(intens_for_kp[, -1]) 
# Creating QFeatures data
cptac <- readQFeatures(intens_for_kp, ecol = 2:10, name = "proteins", fnames = "symbol")
cptac$group <- c(rep("AB16_METTL1_KO", 3),
                 rep("AE4_METTL1_KO", 3), rep("NC2_WT", 3))
colData(cptac)
rowDataNames(cptac)
head(assay(cptac))

# Treating zeros as NA
cptac <- zeroIsNA(cptac, i = seq_along(cptac))
nNA(cptac, i = seq_along(cptac))
#barplot(nNA(cptac, i = seq_along(cptac))$nNAcols$pNA)
table(nNA(cptac, i = seq_along(cptac))$nNArows$nNA)
# Removing rows that have >=7 NAs out of 9
cptac <- filterNA(cptac, i = seq_along(cptac), pNA = 7/9) # 3374 prots left


# IMPUTATION -------------------------------------------------------------------

# Below we are using different imputation methods for the zero imputation in our
# data. Later we will select one method to proceed and will compare it to the 
# results we would get without imputation.

# We have 700+ proteins with NAs overall. 
# Performing imputation and plotting the density of the imputed data. 
cls <- c("black", "red", "blue", "steelblue", "orange")
plot(density(na.omit(assay(cptac))), col = cls[1],  xlim = c(0, 4e+09), ylim = c(0, 5.0e-09))
lines(density(assay(impute(cptac, method = "knn"))), col = cls[2])
lines(density(assay(impute(cptac, method = "zero"))), col = cls[3])
lines(density(assay(impute(cptac, method = "MinDet"))), col = cls[4])
lines(density(assay(impute(cptac, method = "bpca"))), col = cls[5])
legend("topright", legend = c("orig", "knn", "zero", "MinDet", "bpca"),
       col = cls, lwd = 2, bty = "n")

cls <- c("black", "red")
plot(density(na.omit(assay(cptac))), col = cls[1],  xlim = c(0, 2e+09), ylim = c(0, 4.0e-09))
lines(density(assay(impute(cptac, method = "knn"))), col = cls[2])
legend("topright", legend = c("orig", "knn"),
       col = cls, lwd = 2, bty = "n")

# MinDet looks the best but gives terrible distribution after log normalization
# with a second peak in lower tail. Same happens for min imputation. However, 
# knn and bpca imputaions do not shift the distribution!

# So for the further analysis we will use BPCA imputation. 

cptac_imp <- impute(cptac, method = "bpca")
cptac <- impute(cptac, method = "zero")

# Counting unique features
cptac <- countUniqueFeatures(cptac, i = "proteins", colDataName = "proteins_counts")


# TRANSFORMATIONS --------------------------------------------------------------

# Log-transformation
cptac <- addAssay(cptac,
                  logTransform(cptac[[1]]),
                  name = "prots_log")
cptac_imp <- addAssay(cptac_imp,
                      logTransform(cptac_imp[[1]]),
                      name = "prots_log")

# Normalization
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


# DE ANALYSIS ------------------------------------------------------------------

# Below we will be using limma to calculate p-values and logFC.
# The limma will output NA in p-values if protein has NA in more than two groups. 

prots <- getWithColData(cptac, "prots_norm")
prots_imputed <- getWithColData(cptac_imp, "prots_norm")

groups <- prots$group
design <- model.matrix(~ 0 + groups)
contrast_matrix <- makeContrasts(
  AB16_METTL1_KO_vs_NC2_WT = groupsAB16_METTL1_KO - groupsNC2_WT,
  AE4_METTL1_KO_vs_NC2_WT = groupsAE4_METTL1_KO - groupsNC2_WT,
  levels = colnames(design)
)

# Zeros (not-imputed)
fit_prot <- lmFit(assay(prots), design)
contrast_fit_prot <- contrasts.fit(fit_prot, contrast_matrix)
eb_fit_prot <- eBayes(contrast_fit_prot)
# Imputed
fit_prot_imp <- lmFit(assay(prots_imputed), design)
contrast_fit_prot_imp <- contrasts.fit(fit_prot_imp, contrast_matrix)
eb_fit_prot_imp <- eBayes(contrast_fit_prot_imp)

# Zeros
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
# Imputed
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


# SAVING THE RESULTS -----------------------------------------------------------

# Saving the full data
write.csv(res_AE4, file = "./reps3nc2_analysis/full_data/res_AE4.csv", row.names = FALSE)
write.csv(res_AE4_imp, file = "./reps3nc2_analysis/full_data/res_AE4_imp.csv", row.names = FALSE)
write.csv(res_AB16_imp, file = "./reps3nc2_analysis/full_data/res_AB16_imp.csv", row.names = FALSE)
write.csv(res_AB16, file = "./reps3nc2_analysis/full_data/res_AB16.csv", row.names = FALSE)

# Filtering the data
res_AE4_imp <- res_AE4_imp %>% filter(adj.P.Val < 0.05 & abs(logFC)>0.5)
res_AE4 <- res_AE4 %>% filter(adj.P.Val < 0.05 & abs(logFC)>0.5)
res_AB16_imp <- res_AB16_imp %>% filter(adj.P.Val < 0.05 & abs(logFC)>0.5)
res_AB16 <- res_AB16 %>% filter(adj.P.Val < 0.05 & abs(logFC)>0.5)

# Saving the filtered data
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

# Samples intersection
inters_AB16_AE4_zeros_names <- intersect(res_AB16$protein, res_AE4$protein)
inters_AB16_AE4_imp_names <- intersect(res_AB16_imp$protein, res_AE4_imp$protein)

# Adding col with mean logFC
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

# Saving the results
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

# Plotting, venn diagrams for intersections

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

