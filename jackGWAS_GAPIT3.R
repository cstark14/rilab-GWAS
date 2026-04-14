devtools::install_github("jiabowang/GAPIT3",force=TRUE)
library(GAPIT)
library(data.table)

# Load numeric genotype
genotype <- read.table("DH_phg_Jack_filtered_numeric.raw", header = TRUE)

# Drop PLINK metadata columns (FID, PAT, MAT, SEX, PHENOTYPE)
genotype_clean <- genotype[, -c(2:6), with = FALSE]
# Column 1 = IID, rest = SNPs

bim <- read.table("DH_phg_Jack_bed.bim", header = FALSE)

snp_info <- data.frame(
  SNP = bim$V2,   # SNP ID
  Chr = bim$V1,   # chromosome
  Pos = bim$V4    # position
)

myGAPIT <- GAPIT(
  Y = blup_df,
  GD = genotype_clean,
  GM = snp_info,
  PCA.total = 5,
  model = c("BLINK", "FarmCPU", "MLM"),
  Multiple_analysis = TRUE
)
