
## Welcome to the 2022 Consortium SNP Array Data Analysis Workshop -- session for population genetic analyses  

# This repository includes scripts and data associated with the practical sessions of population genetic analyses
# Prepared by Honggang Zhao in Sep 2022, with helps from Matt Hare lab. 

# R packages --------------------------------------------------------------

# Install the R package, please skip this block if you have already installed the package in your system

# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("LEA")
# BiocManager::install("vcfR")
# 
# install.packages("devtools")
# #install hierfstat from devtools:
# library(devtools)
# install_github("jgx65/hierfstat")
# 
# install.packages("SNPfiltR")
# install.packages("pheatmap")
# install.packages("ggplot2")
# install.packages("bigsnpr")

# load packages (please uncomments the codes above when any specific package is not installed yet)
library(LEA)
library(vcfR)
library(hierfstat)
library(SNPfiltR)
library(pheatmap)
library(ggplot2)
library(bigsnpr)

## Part0: set up the working directory for RStudio -------------------------

#set up the working directory
# for Windows users, please use Session-Set Working Directory-Choose Directory to find where the folder HSRL_workshop is located and hit open
# See https://github.com/hzz0024/HSRL_workshop for figure illustration
setwd("~/Desktop/HSRL_workshop")

## Part1: Handling SNP array data: VCF filtering and formatting
# 
# In this demo we will use a vcf file produced from SNP array and converted by Axiom Analysis Suite. 
# The file includes 65,893 SNPs (i.e., 66K array markers) and 125 individuals from four populations. 
  
  # Detailed information for four populations,
  # 
  # |     Population                             | Abbr | N  | Salinity    | Origin          |
  # |--------------------------------------------|------|----|-------------|-----------------|
  # | Lloyd harbor, NY, Long Island Sound wild 1 | LIW1 | 31 | High        | Wild population |
  # | Niantic Bay, CT, Long Island Sound wild 2  | LIW2 | 30 | High        | Wild population |
  # | Rutgers NEH 19N1357                        | NEH1 | 32 | High, 18-23 | Selected line   |
  # | Rutgers NEH 20N1                           | NEH2 | 32 | High, 18-23 | Selected line   |

### Step 1: load vcf using SNPfiltR and snpR ----------------------------------------------------------------------------------------------------------------------------------------------

# load packages
library(SNPfiltR)
library(vcfR)

#read in vcf as vcfR
vcfR <- read.vcfR("./Example_data/example_66k_n125.recode.vcf")

### Step 2: Quality filtering for missing data and minor allele counts ----------------------------------------------------------------------------------------------------------------------

# Checking missing data by SNP and the effect of various cutoffs on the missingness of each sample
missing_by_snp(vcfR)

#choose a value that retains an acceptable amount of missing data in each SNP, here we require that each SNP with <5% missing data
vcfR_missing<-missing_by_snp(vcfR, cutoff = .95)
# cutoff is specified, filtered vcfR object will be returned
# 1.66% of SNPs fell below a completeness cutoff of 0.95 and were removed from the VCF

# Write out vcf files for downstream analyses.
vcfR::write.vcf(vcf_missing_mac, "./example_66k_n125_missing95_mac6.vcf.gz")

## Part2: Principal component analysis (PCA) 
library(SNPfiltR)
library(vcfR)

#read in vcf as vcfR
vcfR <- read.vcfR("./Example_data/example_66k_n125_missing95_mac6_hwe_LD_clump.recode.vcf")

#generate popmap file. Two column popmap with 'id' and 'pop'
popmap<-data.frame(id=colnames(vcfR@gt)[2:length(colnames(vcfR@gt))],pop=substr(colnames(vcfR@gt)[2:length(colnames(vcfR@gt))], 1,4))

# Check the first few content in popmap
head(popmap)
# Perform PCA analyasis
assess_missing_data_pca(vcfR=vcfR, popmap = popmap, thresholds=NULL,clustering = FALSE)

## Part3: Admixture analysis -----------------------------------------------

# lead LEA package
library(LEA)

#change vcf to geno 
LEA::vcf2geno("./Example_data/example_66k_n125_missing95_mac6_hwe_LD_clump.recode.vcf",
              output.file = "example_66k_n125_missing95_mac6_hwe_LD_clump.geno")
# - number of detected individuals: 125
# - number of detected loci:    1000

# modeling ancestry proportions for different K: from K=1 to K=10
obj <- snmf("example_66k_n125_missing95_mac6_hwe_LD_clump.geno", K = 1:10, ploidy = 2,
            entropy = T, CPU =4, project = "new")

# Find the best K from cross-entropy
plot(obj, col = "blue4", cex = 1.4, pch = 19) #---best is 6 here

#choose the best LEA run
best = which.min(cross.entropy(obj, K = 2))

# Plot ancestry proportions across samples
barchart(obj, K=2,run=best,border=T,space=0,
         col=c("grey","orange"), lab=tab_pop$pop,
         xlab = "Individuals", ylab = "Ancestry proportions (K=2)") -> bp

axis(1, at = 1:length(bp$order), 
     labels = popmap[bp$order, "id"], las = 3, 
     cex.axis = .6)

## Part4: Fst statistics ---------------------------------------------------
# load packages
library(hierfstat)
library(pheatmap)
library(vcfR)

# load vcf file and convert it to genind format
vcf_file = "./Example_data/example_66k_n125_missing95_mac6_hwe_LD_clump.recode.vcf"
vcf <- read.vcfR(vcf_file, verbose = FALSE)
df <- vcfR2genind(vcf)
df@pop <- factor(popmap$pop)
# calculate pairwise FST using Weir and Cockerham (1984)
pairwise_fst <- genet.dist(df, method = "WC84") # Estimates pairwise FSTs according to Weir and Cockerham (1984)
# convert the output into matrix
plot_dt <- as.matrix(pairwise_fst)
# check Fst matrix
plot_dt

# plot the Fst matrix with heatmap
pheatmap(plot_dt, display_numbers = T, cellwidth=50, cellheight=40, main="Pairwise FST")

## Part5: Genetic diversity (Heterozygosity and Allelic Richness) ----------

# Calculate the allelic richness
Arich <- allelic.richness(df,min.n=NULL,diploid=TRUE)
ind_mean <- colMeans(x=Arich$Ar, na.rm = TRUE)
ind_mean
#   LIW1     LIW2     NEH1     NEH2 
# 1.981497 1.981681 1.961553 1.929469

# Calculate the obverved and expected heterozygosity
basicstat <- basic.stats(df, diploid = TRUE, digits = 2) 
names(basicstat)

# Obverved heterozygosity
Ho <- colMeans(x=basicstat$Ho, na.rm = TRUE)
Ho
#   LIW1    LIW2    NEH1    NEH2 
# 0.28779 0.29281 0.29348 0.29587 

# Expected heterozygosity
He <- colMeans(x=basicstat$Hs, na.rm = TRUE)
He
#   LIW1    LIW2    NEH1    NEH2 
# 0.31269 0.31558 0.30919 0.30327 







