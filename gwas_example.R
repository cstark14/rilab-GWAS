
library(bigsnpr)


# ------------------------------------------------------------------------\
# 1) Quality control vcf --------
# ------------------------------------------------------------------------\

# Uses bigsnpr to quality control VCF (or other type of file)
# minor allele frequency (MAF) especially important to control for GWAS
# This is just a wrapper around plink so can just plink on its on if that's easier

bedFile <- snp_plinkQC(plink.path = "~/bin/plink2",
                       file.type = "--vcf",
                       prefix.in = "/scratch/geno/all_ZeaSyn_snps_10M",
                       prefix.out = "/scratch/geno/out/1.INTERMEDIATE_10M_all_ZeaSyn_snps_genomissing.1",
                       hwe = 0,
                       mind = 1,
                       maf = 0, # minimum maf
                       geno = .1, # maximum missing 10%
                       extra.options = "--not-chr 0 --allow-extra-chr --min-alleles 2 --max-alleles 2")



bedFile <- snp_plinkQC(plink.path = "~/bin/plink2",
                       file.type = "--bfile",
                       prefix.in = "/scratch/geno/out/1.INTERMEDIATE_10M_all_ZeaSyn_snps_genomissing.1",
                       prefix.out = "/scratch/geno/out/1.10M_ZeaSyn_filtered_maf0.05.maxMissing0.1_only.PGRPphenotyped.lines",
                       hwe = 0,
                       mind = .1,
                       maf = .05, # minimum maf
                       geno = .1, # maximum missing 10%
                       extra.options = "--const-fid 0 --allow-extra-chr --keep /scratch/geno/Z.3.Zeasynthetic_lines_phenotyped.txt --max-maf .95")

bigSNPfile <- snp_readBed(bedfile = bedFile, 
                          backingfile = "/scratch/geno/out/1.10M_ZeaSyn_filtered_maf0.05.maxMissing0.1_only.PGRPphenotyped.lines")


# Will need to deal with NA's in some way
# This impute NA's to major allele
snp <- snp_attach("/scratch/geno/out/1.10M_ZeaSyn_filtered_maf0.05.maxMissing0.1_only.PGRPphenotyped.lines.rds")

# is it 0,1,2 coded? 
# t <- big_counts(snp$genotypes)
# t[,1]

snp$genotypes <- snp_fastImputeSimple(snp$genotypes, method = "mode")

snp_subset(snp,
           backingfile = "/scratch/geno/out/1.10M_ZeaSyn_filtered_maf0.05.maxMissing0.1_only.PGRPphenotyped.lines_NAimp2maj")


# write out bed/bim/fam files
snp <- snp_attach("/scratch/geno/out/1.10M_ZeaSyn_filtered_maf0.05.maxMissing0.1_only.PGRPphenotyped.lines_NAimp2maj.rds")

snp_writeBed(snp, bedfile = "/scratch/geno/out/1.10M_ZeaSyn_filtered_maf0.05.maxMissing0.1_only.PGRPphenotyped.lines_NAimp2maj.bed")


# If bonferroni is too stringent
# get number of independent snps for setting less stringent threshold
obj.bed <- bed("/scratch/gwas_in/1.WiDiv_942g_AGPv4_imputed_maf0.05_maxmaf0.95_masMissing0.1_onlyphenotypedlines.bed")
kept.snp.inds <-
  bed_clumping(obj.bed,
               thr.r2 = .2)

# record this number
effective.bonf <- .05 / length(kept.snp.inds)

# ------------------------------------------------------------------------\
# 2) determine number of PC's --------
# ------------------------------------------------------------------------\

# To control for population structure, need to use first N principal components of the snp matrix in GWAS model
# Below will work if the bigsnp object was created in the previous step. 
# Otherwise, there are other packages that can create similar outputs. 
# Look for elbow of scree plot and consider the percent variance explained by each PC.
# I used the first 4 PC's for my population. 

snpFilePath <- paste0("/scratch/geno/out/1.10M_ZeaSyn_filtered_maf0.05.maxMissing0.1_only.PGRPphenotyped.lines_NAimp2maj.rds")
snp <- snp_attach(snpFilePath)

G <- snp$genotypes
CHR <- snp$map$chromosome
POS <- snp$map$physical.pos

#Change chromsomes to numeric
CHRN <- as.numeric(gsub("Chr_","",CHR))

obj.svd <- bigsnpr::snp_autoSVD(G = G,infos.chr = CHRN,infos.pos = POS,thr.r2 = 0.5)
saveRDS(obj.svd, file = "/scratch/geno/pcs/2.10M_ZeaSyn_10PCs.SVD_object.rds")

PCs <- obj.svd
plot(PCs)
ggsave("/scratch/geno/pcs/2.10M_ZeaSyn_10PCs.scree.png",
       device = "png")

plot(obj.svd, type = "loadings", loadings = 1:10, coeff = 0.4)
ggsave("/scratch/geno/pcs/2.10M_ZeaSyn_10PCs.loadings.png",
       device = "png")

plot(obj.svd, type = "scores")
ggsave("/scratch/geno/pcs/2.10M_10PCs.scores.png", device = "png")

bedfile <- "/scratch/geno/out/1.10M_ZeaSyn_filtered_maf0.05.maxMissing0.1_only.PGRPphenotyped.lines_NAimp2maj.bed"
obj.bed <- bed(bedfile)
obj.svd <- readRDS("/scratch/geno/pcs/2.10M_ZeaSyn_10PCs.SVD_object.rds")
ind.keep <- attr(obj.svd, "subset")

rssq <- bigsnpr:::prod_and_rowSumsSq(
  obj.bed, ind_row = rows_along(obj.bed), ind_col = ind.keep,
  center = obj.svd$center, scale = obj.svd$scale,
  V = matrix(0, length(ind.keep), 0)  # skip product here
)[[2]]
sum(rssq)  # squared frobenius norm of genotype matrix for subset (= total variance explained)
var_exp <- obj.svd$d^2 / sum(rssq)
# signif(var_exp, 2)
# round(cumsum(var_exp), 3)

# plot a scree-like plot of pve
plot.df <- data.frame(pc.num = c(1:10), pve = var_exp * 100)
ggplot(plot.df, aes(x = pc.num, y = pve)) +
  geom_point() +
  geom_line() +
  theme_light() +
  scale_x_continuous(breaks = c(1:10)) +
  labs(x = "Component Number",
       y = "Percent Variance Explained",
       title = "ZeaSyn 600k, Pop Structure Variance Explained by PC's")
ggsave("/scratch/geno/pcs/2.10M_ZeaSyn_percent_variance_explained.png",
       device = "png")

# plot similar but cumulative variance
plot.df <- data.frame(pc.num = c(1:10), pve = round(100 * cumsum(var_exp), 3))
ggplot(plot.df, aes(x = pc.num, y = pve)) +
  geom_point() +
  geom_line() +
  theme_light() +
  scale_x_continuous(breaks = c(1:10)) +
  labs(x = "Component Number",
       y = "Cumulative Percent Variance Explained",
       title = "ZeaSyn 600k, Cumulative Pop Structure Variance Explained by PC's")
ggsave("/scratch/geno/pcs/2.10M_ZeaSyn_cumulative_percent_variance_explained.png",
       device = "png")


# ------------------------------------------------------------------------\
# 3) run gwas --------
# ------------------------------------------------------------------------\

library(rMVP)
library(dplyr)


out.prefix <- "/scratch/gwas_out/2.rMVP_10M_ZeaSyn_filtered_maf0.05.maxMissing0.1_only.PGRPphenotyped.lines_NAimp2maj"

MVP.Data(fileBed = "/scratch/gwas_in/1.10M_ZeaSyn_filtered_maf0.05.maxMissing0.1_only.PGRPphenotyped.lines_NAimp2maj",
         filePhe = "/scratch/gwas_in/Pheno_Subset_to_test_rmvp.csv",
         fileKin = TRUE,
         filePC = TRUE,
         sep.phe = ",",
         priority = "speed",
         out = out.prefix)


out.prefix <- "/scratch/gwas_out/2.rMVP_10M_ZeaSyn_filtered_maf0.05.maxMissing0.1_only.PGRPphenotyped.lines_NAimp2maj"

# inputs
pheno <- read.table(paste0(out.prefix, ".phe"), header = T)
geno <- attach.big.matrix(paste0(out.prefix, ".geno.desc"))
map <- read.table(paste0(out.prefix, ".geno.map"), header = T)
pcs <- attach.big.matrix(paste0(out.prefix, ".pc.desc"))
kin <- attach.big.matrix(paste0(out.prefix, ".kin.desc"))

# Determined in step 2
num.pcs <- 4


# farm cpu specific parameters 
# p.threshold <- NA
# QTN.threshold <- NA
# method.bin <- 'static'
# bin.size <- NA
# bin.selection <- NA
# maxLoop <- NA

# Runs 3 gwas models 
# 1) GLM = generalized linear model
# 2) FarmCPU 
# 3) MLM = mixed linear model 

# Differences in models come down to how the control for population structure. 
# GLM controls the least, MLM the most and FarmCPU somewhere in the middle.
# Follows that GLM will have most false positives/fewest false negatives, MLM the oppposite and FarmCPU somewhere in the middle. 

# FarmCPU is also unique in that it is a multi-locus model so the subsequent manhattan plots will look
# distinct in that they will just have a single snp as significant instead of the classic "skyscrapers" that gives 
# the plot its name. 

mvp.obj <- 
  MVP(phe = pheno[,c(1,4)],
      geno = geno,
      map = map,
      K = kin,
      nPC.GLM = num.pcs,
      nPC.MLM = num.pcs,
      nPC.FarmCPU = num.pcs,
      priority = "speed",
      method = c("GLM", "MLM", "FarmCPU"),
      outpath = "/scratch/gwas_out/mvp_out_test/",
      file.output = c("pmap", "pmap.signal", "log"))

# make phenotype plot
# MVP.Hist(phe=phenotype, file.type="jpg", breakNum=18, dpi=300)

# bonferroni
bonf.cut <- .05 / nrow(geno) 

# Make manhattans of all models
MVP.Report(mvp.obj, plot.type="m", multracks=TRUE, multitraits = F, threshold= bonf.cut, threshold.col= "black", amplify=TRUE,
           bin.size=1e6, chr.den.col=NULL, signal.col=c("red","green"), signal.cex=c(1,1),
           file.type="jpg",dpi=300)

# manhattan of one trait and model
fcpu.res <- read.csv("Days_to_Anthesis_50_AllExps.FarmCPU.csv") %>% 
  select(-REF, -ALT, -Effect, -SE)

MVP.Report(fcpu.res, plot.type="m", col=c("dodgerblue4","deepskyblue"), LOG10=TRUE, ylim=NULL,
           threshold=c(bonf.cut,1e-4), threshold.lty=c(1,2), threshold.lwd=c(1,1), threshold.col=c("black","grey"), 
           amplify=TRUE,chr.den.col=NULL, signal.col=c("red","green"), signal.cex=c(1,1),
           signal.pch=c(19,19),file.type="jpg",memo="",dpi=300)

# qqplots
MVP.Report(imMVP,plot.type="q",col=c("dodgerblue1", "olivedrab3", "darkgoldenrod1"),threshold=1e6,
           signal.pch=19,signal.cex=1.5,signal.col="red",conf.int.col="grey",box=FALSE,multracks=
             TRUE,file.type="jpg",memo="",dpi=300)
