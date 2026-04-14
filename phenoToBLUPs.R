library(tidyverse)
library(dplyr)
library(lme4)
library(readxl)
options(stringsAsFactors=FALSE)
options(scipen=999)
"%notin%" <- Negate("%in%")

DH_Phenos <- read_excel("~/Downloads/DH_phenotype_JLcombined35_cleanup_output_top10pct_large-plots.xlsx", 
                        na = "NA", skip = 1) %>%
  group_by(MaternalParent) %>%
  mutate(Rep=row_number()) %>%
  ungroup() %>%
  rename(`KW100` = `100-KW`) %>%
  mutate(BR=ifelse(BR=='N',0,1)) %>%
  mutate(ME=ifelse(BR=='N',0,1)) %>%
  dplyr::select(-Source)

### lme4 version 
library(lme4)

traits <- colnames(DH_Phenos)[3:(ncol(DH_Phenos)-1)]

# Your data probably looks something like this:
# Line | Env | Rep | Yield | PlantHeight | EarLength
# B73  | Env1|  1  | 12.4  |    145.2    |   18.3

blup_list <- lapply(traits, function(trait) {
  print(trait)
  #formula <- as.formula(paste0(trait, " ~ (1|MaternalParent)"))
  formula <- paste0(trait, " ~ (1|MaternalParent)")
  
  print(formula)
  model <- lmer(formula, data = DH_Phenos, REML = TRUE)
  
  blups <- ranef(model)$MaternalParent
  colnames(blups) <- trait
  blups
})

# Combine
#blup_df <- do.call(cbind, blup_list)
blup_df <- Reduce(function(x, y) merge(x, y, by = "row.names", all = TRUE) %>%
                    column_to_rownames("Row.names"),
                  blup_list)
# Format for PLINK2
pheno_df <- data.frame(
  FID = unique(DH_Phenos$`Sum25Family #`),
  IID = rownames(blup_df),
  blup_df
)

write.table(pheno_df, "~/Desktop/DH_phg_Jack_phenotypes_BLUPs.pheno", 
            row.names = FALSE, quote = FALSE, sep = "\t")

gemma_pheno <- pheno_df %>%
  dplyr::select(-FID, -IID)  # drop the ID columns

write.table(gemma_pheno, "DH_phg_Jack_phenotypes_BLUPs_gemma.pheno",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
writeLines(traits, "trait_names.txt")
