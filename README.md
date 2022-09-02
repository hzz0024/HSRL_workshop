## Welcome to the 2022 Consortium SNP Array Data Analysis Workshop -- session for population genetic analyses  

This repository includes scripts and data associated with the practical sessions of population genetic analyses
Prepared by Honggang Zhao, and Matt Hare in Sep 2022

All the tutorial can be completed with the R/Rstudio. We encourge you to pre-install the R packages that needed for data analyses. These include

[RStudio](https://www.rstudio.com/)

[How to install RStudio on Mac](https://teacherscollege.screenstepslive.com/a/1135059-install-r-and-r-studio-for-mac)

[How to install RStudio on Windows](https://techvidvan.com/tutorials/install-r/#install-r-windows)

R package needed for workshop

[SNPrelate](https://www.bioconductor.org/packages/release/bioc/html/SNPRelate.html)\
[ggplot2](https://ggplot2.tidyverse.org/)

Here is the tentative schedule for the workshop. Within each block we'll keep some flexibility to have questions or discussion.
![schedule](./images_tutorial/schedule.png)


## Part1: Handling SNP array data: filter raw vcf

```R
1. #Replace the original array sample id with read sample id
cd /workdir/hz269/domestication_600K/00_vcf
./1_1_snp_array_format.sh

bcftools query -l genetyped_data_all_samples.vcf > sample_original_name # 842 samples
bcftools reheader -s sample_rename genetyped_data_all_samples.vcf > genetyped_data_all_samples.rename.vcf
mv genetyped_data_all_samples.rename.vcf

 #Replace the chromosome with numbers (1-10)
./1_2_sed_chr.sh

for i in genetyped_data_all_samples.rename.vcf; do
    sed -i.bak 's/NC_035780.1/1/g;s/NC_035781.1/2/g;s/NC_035782.1/3/g;s/NC_035783.1/4/g;s/NC_035784.1/5/g;s/NC_035785.1/6/g;s/NC_035786.1/7/g;s/NC_035787.1/8/g;s/NC_035788.1/9/g;s/NC_035789.1/10/g;s/MT/11/g' $i
done

rm genetyped_data_all_samples.rename.vcf.bak
./1_3_change_ID_array.sh
#!/bin/sh
python3 add_array.py

cat add_array.py
fname = 'genetyped_data_all_samples.rename.vcf'
outname = fname + '.out'

idx = 0
with open(fname, 'r') as f, open(outname, 'w') as w:
    for l in f:
        if l.startswith('#'):
            pass
        else:
            idx += 1
            ss = l.split()
            chrom = ss[0]
            pos = ss[1]
            ID = chrom + '_' + pos
            if ss[2].startswith('AX'):
                ss[2] = ID
                l = '\t'.join(ss)
                l += '\n'
        w.write(l)

# exclude inversions and mtDNA
./1_4_exclude_invers.sh

vcftools --vcf genetyped_data_all_samples.rename.vcf.out --exclude-bed delly.inversions.masked.bed --recode --recode-INFO-all --out genetyped_data_all_samples_nodellyinvers

After filtering, kept 300446 out of a possible 300446 Sites

vcftools --vcf genetyped_data_all_samples_nodellyinvers.recode.vcf --exclude-bed invers.bed --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --recode --recode-INFO-all --out genetyped_data_all_samples_noinvers

After filtering, kept 276327 out of a possible 300446 Sites

# exclude LGF, PCs, VC familes, NYH, and CBW populations, and exclude one individual (DBW1-30) that potentailly mislabeled
1_5_keep_population_n_509.sh

vcftools --vcf genetyped_data_all_samples_noinvers.recode.vcf --keep sample_id_509.txt --recode --recode-INFO-all --out genetyped_data_n_509

VCFtools - 0.1.17
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf genetyped_data_all_samples_noinvers.recode.vcf
	--keep sample_id_509.txt
	--recode-INFO-all
	--out genetyped_data_n_509
	--recode

Keeping individuals in 'keep' list
After filtering, kept 509 out of 842 Individuals
Outputting VCF file...
After filtering, kept 276327 out of a possible 276327 Sites

# check the missing ind
1_6_indmiss.sh

source /programs/miniconda3/bin/activate dDocent-2.8.13
./filter_missing_ind.sh genetyped_data_n_509.recode.vcf genetyped_data_n_509_indmiss

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf genetyped_data_n_509.recode.vcf
	--missing-indv
	--out genetyped_data_n_509_indmiss

After filtering, kept 509 out of 509 Individuals
Outputting Individual Missingness
After filtering, kept 276327 out of a possible 276327 Sites
Run Time = 12.00 seconds



                                          Histogram of % missing data per individual
      400 +---------------------------------------------------------------------------------------------------------+
          |                 *                +                 *                 +                +                 |
          |                 *                                 '*otalmissing' using (bin($1,binwidth)):(1.0) ******* |
      350 |-+               *                                  *                                                  +-|
          |                 *                                  *                                                    |
          |                 *                                  *                                                    |
      300 |-+               *                                  *                                                  +-|
          |                 *                                  *                                                    |
          |                 *                                  *                                                    |
      250 |-+               *                                  *                                                  +-|
          |                 *                                  *                                                    |
      200 |-+               *                                  *                                                  +-|
          |                 *                                  *                                                    |
          |                 *                                  *                                                    |
      150 |-+               *                                  *                                                  +-|
          |                 *                                  *                                                    |
          |                 *                                  *                                                    |
      100 |-+               *                                  *                                                  +-|
          |                 *                                  *                                                    |
          |                 *                                  ************************************                 |
       50 |-+               *                                  *                                  *               +-|
          |******************                                  *                                  *                 |
          |                 *                +                 *                 +                ******************|
        0 +---------------------------------------------------------------------------------------------------------+
        0.025              0.03            0.035              0.04             0.045             0.05             0.055
                                                       % of missing data

The 85% cutoff would be 0.0411867
Would you like to set a different cutoff, yes or no
yes
Please enter new cutoff
0.1
All individuals with more than 10.0% missing data will be removed.

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf genetyped_data_n_509.recode.vcf
	--remove lowDP.indv
	--recode-INFO-all
	--out genetyped_data_n_509_indmiss
	--recode

Excluding individuals in 'exclude' list
After filtering, kept 509 out of 509 Individuals
Outputting VCF file...
After filtering, kept 276327 out of a possible 276327 Sites

# maf and genotype rate filtering
1_7_maf_rate_filter.sh
vcftools --vcf genetyped_data_n_509_indmiss.recode.vcf --maf 0.05 --max-missing 0.95 --recode --recode-INFO-all --out genetyped_data_n_509_maf05_maxmiss095

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
  --vcf genetyped_data_n_509_indmiss.recode.vcf
  --recode-INFO-all
  --maf 0.05
  --max-missing 0.95
  --out genetyped_data_n_509_maf05_maxmiss095
  --recode

After filtering, kept 509 out of 509 Individuals
Outputting VCF file...
After filtering, kept 148063 out of a possible 276327 Sites
Run Time = 40.00 seconds

# call rate missing filtering in each population
1_8_pop_missing.sh

source /programs/miniconda3/bin/activate dDocent-2.8.13
pop_missing_filter.sh genetyped_data_n_509_maf05_maxmiss095.recode.vcf popmap.txt 0.95 17 genetyped_data_n_509_maf05_maxmiss095_popmiss095

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
  --vcf genetyped_data_n_509_maf05_maxmiss095.rename.vcf
  --exclude-positions loci.to.remove
  --recode-INFO-all
  --out genetyped_data_n_509_maf05_maxmiss095_popmiss095
  --recode

After filtering, kept 509 out of 509 Individuals
Outputting VCF file...
After filtering, kept 148063 out of a possible 148063 Sites
Run Time = 30.00 seconds

# filter for HWE
1_9_HWE_by_pop.sh

source /programs/miniconda3/bin/activate dDocent-2.8.13
filter_hwe_by_pop.pl -v genetyped_data_n_509_maf05_maxmiss095_popmiss095.recode.vcf -p popmap.txt -h 0.01 -c 0.5 -o genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe
 
#filter_hwe_by_pop.pl -v <vcffile> -p <popmap> [options]
#
#    Options: -v <vcffile> input vcf file -p <popmap> tab-separated file of
#    samples and population designations -h [hwe] minimum Hardy-Weinberg
#    p-value cutoff for SNPs -c [cutoff] proportion of all populations that a
#    locus can be below HWE cutoff without being filtered -o [out] name of outfile

Processing population: DBW1 (31 inds)
Processing population: DBW2 (32 inds)
Processing population: DBX1 (32 inds)
Processing population: DBX2 (31 inds)
Processing population: DBX3 (31 inds)
Processing population: LIW1 (31 inds)
Processing population: LIW2 (30 inds)
Processing population: MEH2 (32 inds)
Processing population: MEW1 (30 inds)
Processing population: MEW2 (31 inds)
Processing population: NCW1 (32 inds)
Processing population: NCW2 (30 inds)
Processing population: NEH1 (32 inds)
Processing population: NEH2 (32 inds)
Processing population: UMFS (30 inds)
Processing population: UNC1 (20 inds)
Processing population: UNC2 (22 inds)
Outputting results of HWE test for filtered loci to 'filtered.hwe'
Kept 141960 of a possible 148063 loci (filtered 6103 loci)

# access the missing rate, call rate and allele frequency distribution
1_10_summary.sh

vcftools --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf --missing-indv --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe
vcftools --vcf sample_rename_hwe.recode.vcf --missing-site --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe
vcftools --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf --freq2 --max-alleles 2 --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe

########################################
# exclude the union outlier candidates #
########################################

setwd("~/Dropbox/Mac/Documents/HG/Domestication/00_vcf")
system(paste(vcftools," --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf --exclude union_outliers_1575.list --recode --recode-INFO-all --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral", sep=""))
# VCFtools - v0.1.13
# (C) Adam Auton and Anthony Marcketta 2009
# 
# Parameters as interpreted:
#   --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf
# --recode-INFO-all
# --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral
# --recode
# --exclude union_outliers_1575.list
# 
# After filtering, kept 509 out of 509 Individuals
# Outputting VCF file...
# After filtering, kept 140385 out of a possible 141960 Sites
# Run Time = 20.00 seconds

system(paste(plink, " --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral.recode.vcf --allow-extra-chr --make-bed --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral", sep=""))
# PLINK v1.90b6.24 64-bit (6 Jun 2021)           www.cog-genomics.org/plink/1.9/
#   (C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
# Logging to genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral.log.
# Options in effect:
#   --allow-extra-chr
# --make-bed
# --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral
# --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral.recode.vcf
# 
# 16384 MB RAM detected; reserving 8192 MB for main workspace.
# --vcf: 140k variants complete.
# genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral-temporary.bed +
#   genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral-temporary.bim +
#   genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral-temporary.fam
# written.
# 140385 variants loaded from .bim file.
# 509 people (0 males, 0 females, 509 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to
# genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral.nosex .
# Using 1 thread (no multithreaded calculations invoked).
# Before main variant filters, 509 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Total genotyping rate is 0.988849.
# 140385 variants and 509 people pass filters and QC.
# Note: No phenotypes present.
# --make-bed to genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral.bed
# + genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral.bim +
#   genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral.fam ... done.

### produce LD-clumping dataset from neutral data 140,385 SNps

sub_name="genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral"
  
f_bk = paste0(sub_name, ".bk")
if (file.exists(f_bk)) {
  #Delete file if it exists
  file.remove(f_bk)
}

# part 1 SNP clumpping and data preparation
snp_readBed(paste0(sub_name, ".bed"))
# this will create a .rds file
obj.bigSNP <- snp_attach(paste0(sub_name, ".rds"))
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
write.table(keep_snp_ids, file = paste0(sub_name, "_pruned_SNP.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
print(paste0("SNPs after clumpping is: ", length(keep_snp_ids), " out of ", dim(obj.bigSNP$map)[1]))

# generate thinned vcf file
system(paste(vcftools," --vcf ",sub_name,".recode.vcf", " --snps ", sub_name, "_pruned_SNP.txt", " --recode --recode-INFO-all --out ", sub_name, "_pruned", sep=""))

# VCFtools - v0.1.13
# (C) Adam Auton and Anthony Marcketta 2009
# 
# Parameters as interpreted:
#   --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral.recode.vcf
# --recode-INFO-all
# --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral_pruned
# --recode
# --snps genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral_pruned_SNP.txt
# 
# After filtering, kept 509 out of 509 Individuals
# Outputting VCF file...
# After filtering, kept 105660 out of a possible 140385 Sites
# Run Time = 16.00 seconds

system(paste(plink, " --vcf ",sub_name,"_pruned.recode.vcf", " --allow-extra-chr --make-bed --out ", sub_name,"_pruned", sep=""))
# PLINK v1.90b6.24 64-bit (6 Jun 2021)           www.cog-genomics.org/plink/1.9/
#   (C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
# Logging to genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral_pruned.log.
# Options in effect:
#   --allow-extra-chr
# --make-bed
# --out genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral_pruned
# --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral_pruned.recode.vcf
# 
# 16384 MB RAM detected; reserving 8192 MB for main workspace.
# --vcf: 105k variants complete.
# genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral_pruned-temporary.bed
# +
#   genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral_pruned-temporary.bim
# +
#   genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral_pruned-temporary.fam
# written.
# 105660 variants loaded from .bim file.
# 509 people (0 males, 0 females, 509 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to
# genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral_pruned.nosex .
# Using 1 thread (no multithreaded calculations invoked).
# Before main variant filters, 509 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Total genotyping rate is 0.988284.
# 105660 variants and 509 people pass filters and QC.
# Note: No phenotypes present.
# --make-bed to
# genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral_pruned.bed +
#   genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral_pruned.bim +
#   genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_neutral_pruned.fam ...
# done.

### shared outliers among the candidatess
```

## Part1: Handling NGS data: From raw reads to SNPs matrix

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
