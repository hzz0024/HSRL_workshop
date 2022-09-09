## Welcome to the 2022 Consortium SNP Array Data Analysis Workshop -- session for population genetic analyses  

This repository includes scripts and data associated with the practical sessions of population genetic analyses\
Prepared by Honggang Zhao, and Matt Hare in Sep 2022

## Installation of R/RStudio

All the tutorial can be completed with the R/Rstudio. We encourge you to pre-install the R packages and look over related tutorials that needed for data analyses. 

[R](https://www.r-project.org/): R is a free software environment for statistical computing and graphics. It compiles and runs on a wide variety of UNIX platforms, Windows and MacOS. The current R version can be downloaded [here](https://cloud.r-project.org/). Please download the R for your operation system. For MacOS purchased after 2020, please install the R with Apple silicon arm64 build (e.g., [R-4.2.1-arm64.pkg](https://cloud.r-project.org/bin/macosx/big-sur-arm64/base/R-4.2.1-arm64.pkg))

[RStudio](https://www.rstudio.com/): RStudio is an integrated environment for R, a platform helps with code running and data visulization. Note RStudio is not functional without an installation of R. So please install the R before using RStudio. 

[How to install R and RStudio on Mac](https://teacherscollege.screenstepslive.com/a/1135059-install-r-and-r-studio-for-mac)

[How to install R and RStudio on Windows](https://techvidvan.com/tutorials/install-r/#install-r-windows)

## R packages needed for workshop

[SNPRelate](https://www.bioconductor.org/packages/release/bioc/html/SNPRelate.html): SNPrelate is a parallel computing toolset for relatedness and principal component analysis of SNP data. The tutorial of SNPRelate is [here](https://www.bioconductor.org/packages/devel/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html)


How to install:
```sh
# To install the package SNPRelate, you need a current version (>=2.14.0) of R and the R package gdsfmt. After installing R you can run the following commands from the R command shell to install the R package SNPRelate.
# Install the package from Bioconductor repository:
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("gdsfmt")
BiocManager::install("SNPRelate")

# Install the development version from Github. This will install the SNPRelate from the Github. It is the same as "install from sources the packages which need compilation"
library("devtools")
install_github("zhengxwen/gdsfmt")
install_github("zhengxwen/SNPRelate")

# Load the R packages: gdsfmt and SNPRelate
library(gdsfmt)
library(SNPRelate)
```

Tips for SNPrelate installation:
1. Alway type ```library(SNPRelate)``` in the [RStudio Console](https://swcarpentry.github.io/r-novice-inflammation/09-supp-intro-rstudio/index.html) to check if the package has been properly installed or not.
2. Question: what should I do when seeing a message like this: *"Do you want to install from sources the packages which need compilation? (Yes/no/cancel)” and “Update all/some/none? [a/s/n]”*.\
The RStudio asks because the package has updated recently on CRAN but the binary isn't yet available for your OS. Therefore the package manager may need alternative ways to install the SNPrelate. [Yes] should update everything to its latest version (e.g., from the Github), but only if you installed the latest version of pacakge. [No] will ensure that all packages get updated, but not necessarily to their latest versions.  [cancel] will quit the installation process. I would recommand to click "Yes" and check if there is any error message. For question related to *“Update all/some/none? [a/s/n]”*, most the time the ```n``` option works.\
3. Question: what should I do when seeing a message like this: *"compilation failed for package ‘SNPRelate’"*\
You need to check the codes surrounding the error message and figure out where the error comes from. For example, I saw an error message during the SNPRelate installation:
```
ld: warning: directory not found for option '-L/usr/local/gfortran/lib/gcc/11'
ld: library not found for -lquadmath
clang: error: linker command failed with exit code 1 (use -v to see invocation)
make: *** [SNPRelate.so] Error 1
ERROR: compilation failed for package ‘SNPRelate’
```
After some Google search, I found a solution [here](https://github.com/RcppCore/RcppArmadillo/issues/262) and sucessfully install the SNPRelate after reinstall the [gfortran](https://cran.r-project.org/bin/macosx/tools/gfortran-6.1.pkg) for my MacOS (Monterey V12.4, Apple M1).\
4. Always looking for answers by Googling key words in the error message. There is a whole R community in the internet to support you.

***
[LEA](https://bioconductor.org/packages/release/bioc/html/LEA.html): LEA is an R package dedicated to population genomics, landscape genomics and genotype-environment association tests. The tutorial of LEA is [here](http://membres-timc.imag.fr/Olivier.Francois/LEA/files/LEA_github.pdf)

How to install:
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("LEA")

# load LEA
library(LEA)
```

***
[hierfstat](https://cran.r-project.org/web/packages/hierfstat/index.html): The hierfstat package is intended for the analysis of population structure using genetic markers. It is suitable for both haploid and diploid data. In particular, it contains functions to estimate and test hierarchical F-statistics for any number of hierarchical levels.. The tutorial of hierfstat is [here](https://cran.r-project.org/web/packages/hierfstat/vignettes/hierfstat.html)

How to install:
```
You will need the package devtools to be able to install the devel version of hierfstat. To install devtools:

install.packages("devtools")

#install hierfstat from devtools:
library(devtools)
install_github("jgx65/hierfstat")

# load hierfstat
library(hierfstat)
```

Tips for hierfstat installation:
A message will appear during hierfstat installation:

```
Loading required package: usethis
> install_github("jgx65/hierfstat")
Downloading GitHub repo jgx65/hierfstat@HEAD
These packages have more recent versions available.
It is recommended to update all of them.
Which would you like to update?

1: All                           
2: CRAN packages only            
3: None                          
4: httpuv (1.6.5 -> 1.6.6) [CRAN]

Enter one or more numbers, or an empty line to skip updates:
```
An empty line is fine for installation.

***
[vcfR](https://cran.r-project.org/web/packages/vcfR/index.html): vcfR is an R package that facilitates easy manipulation of variant call format (VCF) data. The tutorial of vcfR is [here](https://cran.r-project.org/web/packages/vcfR/vignettes/intro_to_vcfR.html)

How to install:
```
# install the BiocManager from CRAN using the base R install.packages() function
install.packages("BiocManager")

# install the vcfR package from bioconductor using BiocManager::install()
BiocManager::install("vcfR")

# load vcfR
library(vcfR)
```

***
[SNPfiltR](https://devonderaad.github.io/SNPfiltR/index.html): SNPfiltR is an R package to streamline and automate the process of choosing appropriate filtering parameters for next-gen SNP datasets. The tutorial of SNPfiltR is [here](https://devonderaad.github.io/SNPfiltR/articles/reproducible-vignette.html)

How to install:
```
#Install current release from CRAN
install.packages("SNPfiltR")

#Install current development version directly from GitHub
library(devtools)
install_github("DevonDeRaad/SNPfiltR")

# load SNPfiltR
library(SNPfiltR)
```

## R Visualization tools 

***
[ComplexHeatmap](https://www.bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html): Complex heatmaps are efficient to visualize associations between different sources of data sets and reveal potential patterns. Here the ComplexHeatmap package provides a highly flexible way to arrange multiple heatmaps and supports various annotation graphics. The book of ComplexHeatmap is [here](https://jokergoo.github.io/ComplexHeatmap-reference/book/)

How to install:
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")

# load ComplexHeatmap
library(ComplexHeatmap)
```

***
[circlize](https://cran.r-project.org/web/packages/circlize/index.html): circlize package provides an implementation of circular layout generation in R as well as some visuilization options. The book of circlize is [here](https://jokergoo.github.io/circlize_book/book/)

How to install:
```
#The package can be installed from CRAN:
install.packages("circlize")

# directly from GitHub:
devtools::install_github("jokergoo/circlize")

# load circlize
library(circlize)
```

***
[ggplot2](https://ggplot2.tidyverse.org/)

How to install:
```
# The easiest way to get ggplot2 is to install the whole tidyverse:
install.packages("tidyverse")

# Alternatively, install just ggplot2:
install.packages("ggplot2")

# Or the development version from GitHub:
# install.packages("devtools")
devtools::install_github("tidyverse/ggplot2")

# load ggplot2
library(ggplot2)
```

## Optional package for linkage disequilibrium (LD) clumping

[bigsnpr](https://cran.r-project.org/web/packages/bigsnpr/readme/README.html)(optional): bigsnpr is an R package for the analysis of massive SNP arrays, primarily designed for human genetics. The tutorial of bigsnpr is [here](https://privefl.github.io/bigsnpr/index.html)

```
# install.packages("remotes")
remotes::install_github("privefl/bigsnpr")
or for the CRAN version

install.packages("bigsnpr")

# load bigsnpr
library(bigsnpr)
```

***

Within each block below we'll keep some flexibility to have questions or discussion.

## Part0: set up the working directory for RStudio

```
#set up the working directory
setwd("~/Desktop/HSRL_workshop")
```

Or click Session-Set Working Directory-Choose Directory and direct to HSRL_workshop in the Desktop

`Mac`
![result](./Seesion_Mac.png)

`Windows`
![result](./Seesion_Windows.png)

## Part1: Handling SNP array data: VCF filtering and formatting


The Variant Call Format (VCF) file is a data format produced by variant calling software (e.g. Axiom Analysis Suite, GATK, FreeBayes, SAMtools). It contains the information for polymorphic loci (variants) present in the sample or population. The variants can be single nucleotide polymorphism (SNP) or a stretch of insertions or deletions (INDEL). In the VCF file, the variant data is normally represented by 8 columns (#CHROM, POS, ID, REF, ALT, QUAL, FILTER and INFO). The INFO columns contain additional information about the variants, for details of the columns headers please see [here](https://www.reneshbedre.com/blog/vcf-fields.html). 

In this demo we will use a vcf file produced from SNP array and converted by Axiom Analysis Suite. The file includes 65,893 SNPs (i.e., 66K array markers) and 125 individuals. 

Four populations include,

|     Population                             | Abbr | N  | Salinity    | Origin          |
|--------------------------------------------|------|----|-------------|-----------------|
| Lloyd harbor, NY, Long Island Sound wild 1 | LIW1 | 31 | High        | Wild population |
| Niantic Bay, CT, Long Island Sound wild 2  | LIW2 | 30 | High        | Wild population |
| Rutgers NEH 19N1357                        | NEH1 | 32 | High, 18-23 | Selected line   |
| Rutgers NEH 20N1                           | NEH2 | 32 | High, 18-23 | Selected line   |

### Step 1:

#### load vcf using SNPfiltR and snpR

```r
library(SNPfiltR)
#This is SNPfiltR v.1.0.0

#Detailed usage information is available at: devonderaad.github.io/SNPfiltR/ 

#If you use SNPfiltR in your published work, please cite the following papers: 

#DeRaad, D.A. (2022), SNPfiltR: an R package for interactive and reproducible SNP filtering. Molecular Ecology Resources, 00, 1-15. http://doi.org/10.1111/1755-0998.13618 

#Knaus, Brian J., and Niklaus J. Grunwald. 2017. VCFR: a package to manipulate and visualize variant call format data in R. Molecular Ecology Resources, 17.1:44-53. http://doi.org/10.1111/1755-0998.12549
library(vcfR)
#   *****       ***   vcfR   ***       *****
#   This is vcfR 1.13.0 
#     browseVignettes('vcfR') # Documentation
#     citation('vcfR') # Citation
#   *****       *****      *****       *****
```

```r
#set up the working directory

```

```r
#read in vcf as vcfR
vcfR <- read.vcfR("~/example_66k_n125.recode.vcf")
```







```R
vcftools  = "/Users/HG/Dropbox/Mac/Documents/HG/Github/BioinfoTools/vcftools_0.1.13/bin/vcftools";
```

First is to filter vcf to retain SNPs with minor allele frequency of 0.05 and call rate of 0.95.

```R
system(paste(vcftools," --vcf example_66k_n125.recode.vcf --maf 0.05 --max-missing 0.95 --recode --recode-INFO-all --out example_66k_n125_maf05_maxmissing95", sep=""))

# VCFtools - v0.1.13
# (C) Adam Auton and Anthony Marcketta 2009

# Parameters as interpreted:
# 	--vcf example_66k_n219.vcf
# 	--recode-INFO-all
# 	--maf 0.05
# 	--max-missing 0.95
# 	--out example_66k_n219_maf05_maxmissing95
# 	--recode

# After filtering, kept 219 out of 219 Individuals
# Outputting VCF file...
# After filtering, kept 57570 out of a possible 65893 Sites
# Run Time = 5.00 seconds
```

Next we will filter the SNPs based on Hardy-Weinberg equilibrium 

```R
# filter for HWE
system(paste("./filter_hwe_by_pop.pl -v example_66k_n219_maf05_maxmissing95.recode.vcf -p popmap.txt -h 0.01 -c 0.5 -o example_66k_n219_maf05_maxmissing95_hwe"))
 
#filter_hwe_by_pop.pl -v <vcffile> -p <popmap> [options]
#
#    Options: -v <vcffile> input vcf file -p <popmap> tab-separated file of
#    samples and population designations -h [hwe] minimum Hardy-Weinberg
#    p-value cutoff for SNPs -c [cutoff] proportion of all populations that a
#    locus can be below HWE cutoff without being filtered -o [out] name of outfile

# Processing population: LIW1 (31 inds)
# Processing population: LIW2 (30 inds)
# Processing population: NEH1 (32 inds)
# Processing population: NEH2 (32 inds)
# Outputting results of HWE test for filtered loci to 'filtered.hwe'
# Kept 53978 of a possible 55515 loci (filtered 1537 loci)
```
 
Then we evaluate the vcf file for indiviaul call rate, SNP call rate, and the allele frequency distribution

```R
# evaluate the invidual missing rate, SNP call rate and allele frequency distribution
system(paste(vcftools," --vcf example_66k_n219_maf05_maxmissing95_hwe.recode.vcf --missing-indv --out example_66k_n219_maf05_maxmissing95_hwe")
system(paste(vcftools,"  --vcf example_66k_n219_maf05_maxmissing95_hwe.recode.vcf --missing-site --out example_66k_n219_maf05_maxmissing95_hwe")
system(paste(vcftools,"  --vcf example_66k_n219_maf05_maxmissing95_hwe.recode.vcf --freq2 --max-alleles 2 --out example_66k_n219_maf05_maxmissing95_hwe")
```

Usually we want to perform population analyses on indepedent and neutral SNP, here we will perform a LD-clumping step.

```R
# LD clumping 
f_name="example_66k_n219_maf05_maxmissing95_hwe"
f_bk = paste0(f_name, ".bk")
if (file.exists(f_bk)) {
  #Delete file if it exists
  file.remove(f_bk)
}

snp_readBed(paste0(f_name, ".bed"))
# this will create a .rds file
obj.bigSNP <- snp_attach(paste0(f_name, ".rds"))
G <- obj.bigSNP$genotypes
SNPs <- obj.bigSNP$map$marker.ID
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
# check if there is any missing values as NA
#big_counts(G, ind.col = 1:dim(G)[1]) # normally the data include missing values
# genotype imputation
G <- snp_fastImputeSimple(G, method = c("mean0"), ncores = 8) # mean0 is based on rounded mean
#big_counts(G, ind.col = 1:dim(G)[1]) # check if NAs are 0
# LD clumping using r2 = 0.2
newpc <- snp_autoSVD(G, infos.chr = CHR, infos.pos = POS, thr.r2 = 0.2, size = 10) # size is the window size of 10K
# extract SNPs after clumpping
which_pruned = attr(newpc, 'subset')
keep_snp_ids = SNPs[which_pruned]
write.table(keep_snp_ids, file = paste0(f_name, "_thinned_SNP.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
print(paste0("SNPs after clumpping is: ", length(keep_snp_ids), " out of ", dim(obj.bigSNP$map)[1]))

# generate thinned vcf file
system(paste(vcftools," --vcf ",f_name,".recode.vcf", " --snps ", f_name, "_thinned_SNP.txt", " --recode --recode-INFO-all --out ", f_name, "_thinned", sep=""))

# VCFtools - v0.1.13
# (C) Adam Auton and Anthony Marcketta 2009

# Parameters as interpreted:
# 	--vcf example_66k_n219_maf05_maxmissing95_hwe.recode.vcf
# 	--recode-INFO-all
# 	--out example_66k_n219_maf05_maxmissing95_hwe_thinned
# 	--recode
# 	--snps example_66k_n219_maf05_maxmissing95_hwe_thinned_SNP.txt

# After filtering, kept 219 out of 219 Individuals
# Outputting VCF file...
# After filtering, kept 46008 out of a possible 55462 Sites
# Run Time = 3.00 seconds

```
 
## Part2: Principal component analysis (PCA) 

```R
library(SNPRelate)
library(ggplot2)
# load vcf file 
vcf.fn <- "example_66k_n219_maf05_maxmissing95_hwe_thinned.recode.vcf"
# VCF => GDS
snpgdsVCF2GDS(vcf.fn, "example_66k_n219_maf05_maxmissing95_hwe_thinned.recode.gds", method="biallelic.only")
# summary
snpgdsSummary("example_66k_n219_maf05_maxmissing95_hwe_thinned.recode.gds")
# Open the GDS file
genofile <- snpgdsOpen("example_66k_n219_maf05_maxmissing95_hwe_thinned.recode.gds")
# PCA function
pca <- snpgdsPCA(genofile,autosome.only=FALSE)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
# For data sets with a handful of dimensions, one typically retains the first few PCs. 
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],    # the second eigenvector
                  EV4 = pca$eigenvect[,4],    # the second eigenvector
                  stringsAsFactors = FALSE)
print(tab)
# output the tab contents for modification
write.table(tab, "sample_eigen_219.txt", row.names=F, sep="\t", quote=F,col.names=T)
tab_pop = read.delim("sample_eigen_219.txt", header = TRUE, sep='\t')
tab_pop$pop = stringr::str_remove(tab_pop$sample.id, "-[0-9]+")

ggplot(tab_pop, aes(x = EV1, y = EV2)) + 
  geom_point(size=2, aes(color=pop), alpha = 0.8)+
  guides(fill = guide_legend(override.aes=list(shape=17)))+
  scale_x_continuous(paste("PC1 (",round(pc.percent[1],3),"%", ")",sep="")) + 
  scale_y_continuous(paste("PC2 (",round(pc.percent[2],3),"%",")",sep=""))
```

![result](./p2.jpeg)

## Part3: Population structure

```R

library(LEA)

#change vcf to geno 
LEA::vcf2geno("example_66k_n219_maf05_maxmissing95_hwe_thinned.recode.vcf",
              output.file = "example_66k_n219_maf05_maxmissing95_hwe_thinned.geno")

#modeling ancestry proportions for different K: from K=1 to K=10
obj <- snmf("example_66k_n219_maf05_maxmissing95_hwe_thinned.geno", K = 1:10, ploidy = 2,
            entropy = T, CPU =4, project = "new")

# Find the best K from cross-entropy
plot(obj, col = "blue4", cex = 1.4, pch = 19) #---best is 6 here

#choose the best LEA run
best = which.min(cross.entropy(obj, K = 2))

# Plot ancestry proportions across samples
barchart(obj, K=2,run=best,border=NA,space=0,
         col=c("red","yellow","blue"),
         xlab = "Individuals", ylab = "Ancestry proportions (K=2)")
```

![result](./p3.jpeg)

## Part4: Fst statistics

```R
library(hierfstat)
library(vcfR)
library(ComplexHeatmap)
library(circlize)

# load the population information
tab_pop <- read.delim("sample_eigen_219.txt", header = TRUE, sep='\t')
tab_pop$pop = stringr::str_remove(tab_pop$sample.id, "-[0-9]+")
# load vcf file and convert it to genind format
vcf_file = "example_66k_n219_maf05_maxmissing95_hwe_thinned.recode.vcf"
vcf <- read.vcfR(vcf_file, verbose = FALSE)
df <- vcfR2genind(vcf)
df@pop <- factor(tab_pop$pop)
pairwise_fst <- genet.dist(df, method = "WC84") # Estimates pairwise FSTs according to Weir and Cockerham (1984)

plot_dt <- as.matrix(pairwise_fst)

Heatmap(plot_dt, name = "Pairwise Fst",
        col = colorRamp2(c(0, 0.05, 0.1), c("white", "skyblue", "orange")),
        rect_gp = gpar(col = "white", lwd = 2),
        cluster_rows = FALSE, 
        cluster_columns = FALSE)

```

![result](./p4.jpeg)

## cheatsheet for R code
[cheatsheet](https://www.rstudio.com/resources/cheatsheets/)
