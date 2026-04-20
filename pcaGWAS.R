library(ggplot2)
library(data.table)
library(tidyverse)
library(dplyr)
library(qqman)
library(CMplot)

# Load eigenvec and eigenval
eigenvec <- read.table("DH_phg_Jack_pca.eigenvec", header = FALSE)
eigenval <- read.table("DH_phg_Jack_pca.eigenval", header = FALSE)

# Calculate % variance explained
eigenval$pct <- eigenval$V1 / sum(eigenval$V1) * 100
eigenval$PC <- paste0("PC", 1:nrow(eigenval))

colnames(eigenvec) <- c("IID", paste0("PC", 1:(ncol(eigenvec) - 2)))

ggplot(eigenval, aes(x = reorder(PC, -pct), y = pct)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_line(aes(group = 1)) +
  geom_point() +
  labs(x = "Principal Component",
       y = "% Variance Explained",
       title = "PCA Scree Plot") +
  theme_bw()

ggplot(eigenvec, aes(x = PC1, y = PC2, label = IID)) +
  geom_point(alpha = 0.7, size = 2, color = "steelblue") +
  labs(x = paste0("PC1 (", round(eigenval$pct[1], 1), "%)"),
       y = paste0("PC2 (", round(eigenval$pct[2], 1), "%)"),
       title = "PCA Plot") +
  theme_bw()

# # Manhattan plot
# manhattan(gwas_plot,
#           main = "Yield GWAS",
#           suggestiveline = FALSE,
#           genomewideline = FALSE,
#           col = c("steelblue", "orange"))
# abline(h=-log10(0.05/10255162))
# 
# # QQ plot
# qq(gwas_plot$P, main = "Yield Q-Q Plot")
assocFiles <- list.files("outputGemma", pattern = "assoc.txt", full.names = TRUE)
traits <- readLines("trait_names.txt")
if(length(assocFiles)!= length(traits)){
  print("ERROR THE TRAIT LIST DOESN'T MATCH NUMBER OF GEMMA RUNS")
  break
}
for(i in 1:length(traits)){
  currentTrait <- traits[i]
  currentFileBasedOnTrait <- assocFiles[grepl(currentTrait, assocFiles, fixed = TRUE)]
  traitAssoc_gemma <- read.table(file=currentFileBasedOnTrait,header=T)
  
  gwasTable <- data.frame(
    SNP = traitAssoc_gemma$rs,
    CHR = as.numeric(traitAssoc_gemma$chr),
    BP  = traitAssoc_gemma$ps,
    P   = traitAssoc_gemma$p_wald
  ) %>%
    filter(CHR %in% 1:10) %>%
    mutate(FDR = p.adjust(P, method = "BH"))
  
  bonferroniThreshold <- 0.05/nrow(gwasTable)
  
  fdr_threshold <- gwasTable %>%
    filter(FDR < 0.05) 
  if(nrow(fdr_threshold >0)){
    fdr_threshold <- fdr_threshold %>%
      summarise(threshold = max(P)) %>%
      pull(threshold)
    threshold <- c(fdr_threshold,bonferroniThreshold)
    thresholdColor <- c("blue", "red")
  } else{
    threshold <- c(bonferroniThreshold)
    thresholdColor <- c("red")
    print("NO SNPS WITH FDR LESS THAN 0.05")
  }
  
  gwas_plot <- gwasTable %>% dplyr::select(- FDR)
  
  ### CM Plot
  CMplot(gwas_plot,
         plot.type = c("m"),       # manhattan ,"c" + QQ
         #threshold = c(1e-5, 1e-6),
         #threshold.col = c("blue", "red"),
         #threshold=0.05/10255162,
         threshold=threshold,
         threshold.col=thresholdColor,
         amplify = TRUE,
         file = "png",
         main=currentTrait,
         dpi = 300,
         file.output = TRUE,
         file.name=currentTrait)
}

