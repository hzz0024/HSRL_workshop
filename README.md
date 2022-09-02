## Welcome to the 2022 Consortium SNP Array Data Analysis Workshop -- session for population genetic analyses  

This repository includes scripts and data associated with the practical sessions of population genetic analyses
Prepared by Honggang Zhao, and Matt Hare in Sep 2022

All the tutorial can be completed with the R/Rstudio. We encourge you to pre-install the R packages that needed for data analyses. These include

[RStudio](https://www.rstudio.com/)

[How to install RStudio on Mac](https://teacherscollege.screenstepslive.com/a/1135059-install-r-and-r-studio-for-mac)

[How to install RStudio on Windows](https://techvidvan.com/tutorials/install-r/#install-r-windows)

R package needed for workshop

[SNPrelate](https://www.bioconductor.org/packages/release/bioc/html/SNPRelate.html)\n
[ggplot2](https://ggplot2.tidyverse.org/)

Here is the tentative schedule for the workshop. Within each block we'll keep some flexibility to have questions or discussion.
![schedule](./images_tutorial/schedule.png)


## Open Rstudio and load data for data filtering

[AWS for Mac OS X and Linux users](https://github.com/clairemerot/physalia_adaptation_course/blob/2022/AWS_mac.md)

For windows users, please use MobaXterm adn FilZilla following the document sent by Carlo.
https://github.com/clairemerot/physalia_adaptation_course/blob/2021/Connection%20to%20the%20Amazon%20EC2%20service_.pdf

This is a guide for windows but it uses different softwares (putty/winSCP) a bit more complex
[AWS for windows users](https://github.com/clairemerot/physalia_adaptation_course/blob/2022/AWS_windows.md)


## Day1: Handling NGS data: From raw reads to SNPs matrix

Capelin: Data taken from Cayuela et al,2020. Molecular Ecology https://doi.org/10.1111/mec.15499
 
Genome assembly: For this course, we made a dummy assembly of about 90 MB (instead of about 500 MB) and 5 chromosomes (instead of 24).

Raw reads: Obtained by GBS methods (= Rad-seq), sequenced with IonTorrent.

1-1: Getting familiar with Unix environment

1-2: From raw sequences to mapped reads

1-3: Calling variants with Stacks
 
## Day2: Population structure and confounding factors

2-1: Fst statistics with vcftools (optionnal: with Stacks, optional: Pairwise-Fst and isolation-by-Distance )

2-2: Principal component analysis (PCA)

2-3: Population clustering with FastStructure

2-4: Discriminant Analysis of Principal Components (DAPC)

## Day3: Outlier detection - Environmental associations

Data: We focus on 12 population from Canada for which there is almost no geographic structure but environmental variability

3-1: Genetic structure and LD-pruned data

3-2: Outlier of differentiation with two methods (Outflank & Baypass)

3-3: Genotype-Environnement Associations with two methods (Baypass & redundancy analysis)

## Day4: Accounting for Structural Variants

We focus on 12 population from Canada. We recommend that you pick one of the two tutorials (haploblocks by local PCA or CNVs on RAD-seq data)

4-1: Investigating haplotypes blocks ( ~inversions?)

This tutorial include working on local PCA, but also calculation of LD, Fst and observed fraction of heterozygotes which may be useful in other contexts

4-2: Filtering duplicated loci in RAD-seq data ( ~ Copy number variants)

This tutorial show how to filter RAD loci to exclude duplicated ones (keep a reliable dataset for SNP analysis), and then how to analyse the duplicated loci for environmental associations.

4-3: Detecting SV with Delly ??

## Day5: Functional approaches

5-1: SNPeff annotation of SNPs for coding & regulatory regions

5-2: Intersection between SNPs and genes with bedtools

5-3: Gene ontology enrichment

5-4: (Optional) Intersection between CNVs and repeats/TE

## cheatsheet for terminal command lines
![cheatsheet](./images_tutorial/bash_cheatsheet.png)
