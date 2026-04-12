library("sommer")

###################
#Read in & explore data
###################

data(DT_cpdata)

#phenotype data for 363 accessions
DT_cpdata[1:10,1:10] 
dim(DT_cpdata)

#genotype matrix, 2889 SNPs.
GT_cpdata[1:10,1:10] 
dim(GT_cpdata)

#SNP map information, for Manhattan plot.
MP_cpdata[1:10,]    

###################
#Construct G relationship matrix
###################

#Construct the Genomic relationship matrix (similar to rrBLUP)
G_matrix = A.mat(GT_cpdata) # G matrix (genomic relationship matrix)
dim(G_matrix) #equal to number of accessions x number of accessions

###################
#Conduct Principal Components Analysis or Eigenvalue decomposition
###################

#Principal components analysis on marker matrix
pca.results = prcomp(GT_cpdata)
pca <- pca.results 
plot(pca$x[,1], pca$x[,2]) #First two PCs
vars <- apply(pca$x, 2, var)  
props <- vars / sum(vars)
props[1:5]
plot(props, ylab="Fraction Explained", xlab = "Number of Principal Components")

#Comparable to doing eigendecomposition on the G relatiopnship matrix
#(This can be less computationally intensive)
eig.results = eigen(G_matrix)
lambda = eig.results$values
props2 = lambda/sum(lambda)  #first 3 components explain 11% total variance
props2[1:5]
plot(props2, ylab="Fraction Explained", xlab = "Number of Principal Components")

###################
#Conduct a Q + K GWAS
###################

#  General model:
#  y=mu+X*beta+PC+id
#  mu : overall mean
#  X*beta : SNP fixed effect
#  PC : PCs included for Q-type GWASes
#  id : kinship/relationship matrix for K-type GWASes

#Set up sommer gwas model
mix1 = sommer::GWAS(color~1,                    # y~fixed
            random=~vs(id,Gu=G_matrix), # random effect with known covariance matrix G (= genomic relationship matrix) for K-type GWAS
            rcov=~units,                # specify the name of the error term
            data=DT_cpdata,             # where the data resides
            n.PC=3,                     # number of principal components to include as fixed effects (for Q-type GWAS)
            M=GT_cpdata,     #SNP matrix 
            gTerm = "u:id")  #a character indicating the random effect linked to SNP matrix


###################
#Construct a Manhattan Plot (using manhattan())
###################

ms = as.data.frame(mix1$scores)
ms$Locus = rownames(ms)
MP2 = merge(MP_cpdata,ms,by="Locus",all.x = TRUE);
manhattan(MP2, 
          pch=20,
          cex=1, 
          PVCN = "color", 
          fdr.level = 0.05) 
# pch: point shape.
# cex: number indicating the amount by which plotting text and symbols should be scaled relative to the default.
#      1=default, 1.5 is 50% larger.
# PVCN: the "color" column is used for y-axis
# fdr.level: desired FDR threshold

#Add a bonferroni-corrected line using abline(h=-log10(desired p-value threshold/# SNPS))
abline(h=-log10(0.05/2889))

###################
#Construct a Manhattan Plot (using ggplot) : more flexible aesthetically
###################

library(dplyr)
#Re-calculate x-positions for manhattan plot so that chromosome positions lie end-to-end
MP2.gg <- MP2 %>% 
  # Compute chromosome size
  group_by(Chrom) %>% 
  summarise(Chrom_length=max(Position)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(Chrom_length)-Chrom_length) %>%
  select(-Chrom_length) %>%
  # Add this info to the initial dataset
  left_join(MP2, ., by=c("Chrom"="Chrom")) %>%
  # Add a cumulative position of each SNP
  arrange(Chrom, Position) %>%
  mutate(Cumulative_Position = Position + tot)

#Prep x-axis for recalculated x-positions
axisdf = MP2.gg %>% 
  group_by(Chrom) %>% 
  summarize(center=(max(Cumulative_Position) + min(Cumulative_Position) ) / 2 )

library(ggplot2);library(ggthemes)
ggplot(MP2.gg, aes(x=Cumulative_Position, y=color)) +
  geom_point(aes(color=as.factor(Chrom)), alpha=0.8, size=1.3) +
  # choose desired colors
  scale_color_manual(values = rep(c("black", "dark grey"), length(unique(MP2.gg$Chrom)) )) +
  # customize the X axis using the dataframe you made above
  scale_x_continuous(label = axisdf$Chrom, 
                     breaks= axisdf$center) +
  xlab("Chromosome")+
  # customize the Y axis:
  ylab(expression(paste("-log"[10]," (",italic("p"),")")))+
  # add bonferroni threshold
  geom_hline(aes(yintercept = -log10(0.05/2889)), 
             linetype="dashed")+
  # Customize the theme to match your desires
  theme_bw() +
  theme(legend.position="none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size=10),
    axis.text.y = element_text(size=10),
    axis.title.x = element_text(size=10),
    axis.title.y = element_text(size=10))

###################
#Calculate Genomic Inflation Factor
###################

#Function: given marker score, returns p-value
PVAL = function(score){
  pval = sort(10^(-score))
  return(pval)
}

#Function: calculates lambda given p-values
lambda = function(pval){
  chisq = qchisq(1-pval,1)
  lambda = median(chisq)/qchisq(0.5,1)
  return(lambda)
}

#Transform -log10(p) to p
pval = PVAL(MP2$color)

#Calculate lambda byb converting p to to chi-square values
lambda(pval)


###################
#Construct QQ Plots (using qqPlot() from GWASTools)
###################
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("GWASTools")

library(GWASTools)
qqPlot(pval)

###################
#Construct QQ Plots (using ggplot2): more flexible for aesthetics
###################
gg_qqplot <- function(ps, ci = 0.95) {
  #how many p-values are there?
  n  <- length(ps)
  #log observed vs expected values
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  #set up labels
  log10Pe <- expression(paste("Expected -log"[10]," (",italic("p"),")"))
  log10Po <-expression(paste("Observed -log"[10]," (",italic("p"),")"))
  #create ggplot object and add the confidence interval (CI) ribbon
  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    #Add Points, expected line, and
    geom_point(aes(expected, observed), shape = 1, size = 3) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po)
}

#You can adjust themes post-function
gg_qqplot(pval)+
  #base_size controls all text size if you don't want to manually control
  theme_bw(base_size = 10) +
  theme(
    axis.ticks = element_line(size = 0.5),
    panel.grid = element_blank())
