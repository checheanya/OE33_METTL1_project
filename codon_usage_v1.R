########################## Codon Usage Analysis ################################

library(biomaRt)
library(Biostrings)
devtools::install_version("dbplyr", version = "2.3.4")
library(dplyr)
library(ggplot2)

# Loading sequences for all transcripts 
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
transcript_ids <- rownames(contrast_list_trs[[1]]$table)
sequences <- getSequence(id = transcript_ids, type = "ensembl_transcript_id",
                         seqType = "coding", mart = ensembl)

# Selecting up-/down-regulated transcripts and non-altered transcripts for 
# each METTL1-KO vs control comparison 
df_to_plot <- data.frame(txid = rep(transcript_ids, 4),
                         contrast = rep(c('AE4_AB1', 'AB16_AB1', 'AE4_AE21', 'AB16_AE21'), each = length(transcript_ids)),
                         stringsAsFactors = FALSE)

# Function to infer direction for each transcript
direction_func <- function(logFC, PValue) {
  if(abs(logFC) < 0.5 | PValue > 0.05) {
    return('none')
  } else if(logFC < 0) {
    return('down')
  } else {
    return('up')
  }
}

# Initialize an empty list to store direction vectors
dir_vec <- vector("character", length = length(transcript_ids) * 4)

# Fill dir_vec with directions for each transcript and contrast
for (i in 1:4) {
  contrast_table <- contrast_list_trs[[i]]$table
  dir_vec[((i-1) * length(transcript_ids) + 1):(i * length(transcript_ids))] <- 
    sapply(1:nrow(contrast_table), function(tx) direction_func(contrast_table[tx, 'logFC'], contrast_table[tx, 'PValue']))
}

# Assign the calculated directions to the expr column
df_to_plot$expr <- dir_vec

# 61 sense codons
sense_codons <- setdiff(names(GENETIC_CODE), c('TAA','TAG','TGA'))
# m7G affected codons
m7Gcodon = c("GCT", "GCG", "GCA", "AGA", "AAC", "TGC", "GGT", "ATT", "AAG", "AAA",
             "ATG", "TTC", "ACA", "TGG", "TAC", "CCT", "CCG", "CCA", "GTT", "GTG",
             "GTA","ATG", "GGG", "ATA", "TCC", "ACT")

# Calculating METTL1-codons percent for each sequence
percent_vec <- c()
for(j in 64150:nrow(sequences)){
  tr <- sequences[j,]$coding
  sum_m7g_codons <- 0
  sum_nonm7g_codons <- 0
  # counting m7g and non-m7g codons
  if (nchar(tr) >= 6) {
    for(k in seq(4, nchar(tr) - 2, 3)){
      codon <- substr(tr, k, k + 2)
      if(codon %in% m7Gcodon){
        sum_m7g_codons <- sum_m7g_codons + 1
      } else {
        sum_nonm7g_codons <- sum_nonm7g_codons + 1
      }
    }
  }
  # calculating percentage of m7g codons 
  if ((sum_m7g_codons + sum_nonm7g_codons) > 0) {
    percent <- sum_m7g_codons / (sum_m7g_codons + sum_nonm7g_codons)
  } else {
    percent <- 0
  }
  percent_vec <- c(percent_vec, percent)
}

df_to_plot <- filter(df_to_plot, txid %in% sequences$ensembl_transcript_id)
df_to_plot$percent <- percent_vec
df_to_plot <- filter(df_to_plot, percent != 0)
df_to_plot <- filter(df_to_plot, percent != 1)
head(df_to_plot)

# Calculating p-vals for non-altered versus up and down regulated transcripts
pairwise_kruskal <- function(data) {
  p_none_up <- kruskal.test(percent ~ expr, data = data %>% filter(expr %in% c('none', 'up')))$p.value
  p_none_down <- kruskal.test(percent ~ expr, data = data %>% filter(expr %in% c('none', 'down')))$p.value
  return(data.frame(p_none_up = p_none_up, p_none_down = p_none_down))
}

kruskal_results <- df_to_plot %>%
  group_by(contrast) %>%
  do(pairwise_kruskal(.)) %>%
  ungroup() %>%
  mutate(label_up = paste0("none vs up: p = ", signif(p_none_up, 3)),
         label_down = paste0("none vs down: p = ", signif(p_none_down, 3)))

# Plotting
ggplot(df_to_plot, aes(x=contrast, y=percent, fill=expr)) + 
  geom_boxplot() +
  geom_jitter(aes(color=expr), width=0.2, size=1) +
  theme_minimal() +
  labs(title="Percent of m7G codons in the up-, down-regulated transcripts and ones with non-altered expression",
       x="Contrast",
       y="Percent of m7G codons")


