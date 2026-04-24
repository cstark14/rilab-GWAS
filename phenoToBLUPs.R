library(tidyverse)
library(dplyr)
library(lme4)
library(readxl)
#options(stringsAsFactors=FALSE)
options(scipen=999)
"%notin%" <- Negate("%in%")

DH_Phenos <- read_excel("~/Desktop/JackGWAS/DH_phenotype_JLcombined35_cleanup_output_top10pct_large-plots.xlsx", 
                        na = "NA", skip = 1) %>%
  mutate(MaternalParent = gsub(" ","",MaternalParent)) %>%
  group_by(MaternalParent) %>%
  mutate(Rep=row_number()) %>%
  ungroup() %>%
  rename(`KW100` = `100-KW`) %>%
  mutate(BR=ifelse(BR=='N',0,1)) %>%
  mutate(ME=ifelse(ME=='N',0,1)) %>%
  dplyr::select(-Source) 

traits <- colnames(DH_Phenos)[3:(ncol(DH_Phenos)-1)]
binary_traits <- c("BR", "ME")
continuous_traits <- traits[traits %notin% binary_traits]

# Your data probably looks something like this:
# Line | Env | Rep | Yield | PlantHeight | EarLength
# B73  | Env1|  1  | 12.4  |    145.2    |   18.3

BLUPFunctions <- function(){
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
}

BLUEfunctions <- function(){
  
  genotype_ids <- unique(DH_Phenos$MaternalParent)
  
  trait_data <- DH_Phenos %>%
    filter(!is.na(.data[[trait]]))
  
  genotype_ids_trait <- unique(trait_data$MaternalParent)
  
  formula <- as.formula(paste(trait, "~ MaternalParent"))
  model <- glm(formula, data = trait_data, family = binomial)
  
  blues <- predict(model, newdata = list(MaternalParent = genotype_ids_trait),
                   type = "response")
  
  data.frame(MaternalParent = genotype_ids_trait, value = blues) %>%
    rename(!!trait := value)
  
  # Continuous traits
  blue_list_continuous <- lapply(continuous_traits, function(trait) {
    trait_data <- DH_Phenos %>%
      filter(!is.na(.data[[trait]]))
    
    genotype_ids_trait <- unique(trait_data$MaternalParent)
    
    formula <- as.formula(paste(trait, "~ MaternalParent"))
    model <- lm(formula, data = trait_data)
    blues <- predict(model, newdata = list(MaternalParent = genotype_ids_trait))  # fix here
    data.frame(MaternalParent = genotype_ids_trait, value = blues) %>%
      rename(!!trait := value) 
  })
  
  # Binary traits
  blue_list_binary <- lapply(binary_traits, function(trait) {
    trait_data <- DH_Phenos %>%
      filter(!is.na(.data[[trait]]))
    
    genotype_ids_trait <- unique(trait_data$MaternalParent)
    
    formula <- as.formula(paste(trait, "~ MaternalParent"))
    model <- lm(formula, data = trait_data)
    blues <- predict(model, newdata = list(MaternalParent = genotype_ids_trait))  # fix here
    data.frame(MaternalParent = genotype_ids_trait, value = blues) %>%
      rename(!!trait := value)
  })
  
  # Combine
  blue_list <- c(blue_list_continuous, blue_list_binary)
  blue_df <- Reduce(function(x, y) full_join(x, y, by = "MaternalParent"), blue_list)
}

# Format for PLINK2
pheno_df <- data.frame(
  FID = unique(DH_Phenos$`Sum25Family #`),
  IID = blue_df$MaternalParent,
  blue_df %>% select(-MaternalParent)
)

write.table(pheno_df, "DH_phg_Jack_phenotypes_BLUEs.pheno", 
            row.names = FALSE, quote = FALSE, sep = "\t")

gemma_pheno <- pheno_df %>%
  dplyr::select(-FID, -IID)  # drop the ID columns

write.table(gemma_pheno, "DH_phg_Jack_phenotypes_BLUEs_gemma.pheno",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
writeLines(traits, "trait_names.txt")
