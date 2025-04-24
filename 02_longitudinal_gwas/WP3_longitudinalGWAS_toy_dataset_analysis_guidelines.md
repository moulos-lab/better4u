---
title: "[BETTER4U - WP3] Analysis plan for longitudinal BMI/weight trajectories (toy
  dataset) v1.0"
author: "Lampros Bouranis, Harokopio University of Athens"
output:
  word_document:
    number_sections: yes
    fig_caption: yes
  html_document:
    df_print: paged
  pdf_document: default
bibliography: WP3_longitudinalGWAS_analysisplan.bib
csl: "american-statistical-association.csl"
documentclass: article
fontsize: 12pt
geometry: margin = 1in
---

<style type="text/css">

body, td {
   font-size: 14px;
}
code.r{
  font-size: 12px;
}
pre {
  font-size: 20px
}
</style>

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

<style>body {text-align: justify}</style>

# Contact details {-}

- Lampros Bouranis (lampros.bouranis@hua.gr)

# Introduction

This document focuses on **Task 3.4** "Mapping of genetic variants associated with fast or slow weight gain trajectories in the life course". 
Following the data quality control (QC) implemented during **Task 3.2** "Trans-ethnic meta-analysis of the joint genetic architecture of BMI and associated traits in diverse populations" (please refer to the document "Analysis Plan for Genetic Data - BETTER4U WP3 - version 1.1"), each individual study will perform longitudinal genome-wide association analyses (GWAS) for body weight and Body Mass Index (BMI) and will provide summary results for meta-analysis. 
Results files will be deposited to a central repository where the meta-analysis will be performed. Based on sample data provided by HUA, the workflow presented below aims at quantifying the effects of genetic and lifestyle (environmental) factors on BMI and weight change over the time.

# Considerations

Due to the diversity in measurements between cohorts, please only include data measured using the following measurement units:

1. For BMI, please use $kg/m^2$ and two decimal places (e.g. $23.87$).

2. For weight, please use as unit $kg$ with two decimal places. Only perform analyses in the sample of age greater than 18 years old.

3. For height, please use as unit $m$ in two decimals.

4. For gender, please use the code: *1 = male*, *2 = female* as per the BETTER4U codebook.

To implement the proposed workflow, please consider the [sample phenotype file](https://drive.google.com/drive/folders/1Y7jNcJmayhF9PdcmvBACFktSfFHcabWA?usp=drive_link).

# Requirements

## General assumptions

This document complements ["BETTER4U Weight change analyses guidelines"](https://docs.google.com/document/d/1DY6mJ53g5pCpgH-hz5FfMqiktL34Cx-S7qXMGwGRqFg/edit?usp=sharing). The proposed workflow is based on the [sample datasets](https://drive.google.com/drive/folders/1Y7jNcJmayhF9PdcmvBACFktSfFHcabWA?usp=sharing) provided by the HUA team. It is assumed that the genotypic data have undergone the quality control (QC) process described in ["Analysis Plan for Genetic Data - BETTER4U WP3 - version 1.1"](https://drive.google.com/file/d/1nvzHtfC7oOUYQYSg2J9M3T2Ek5dHfdtu/view?usp=drive_link), resulting in a single PLINK file with all autosomal chromosomes. For the purpose of this demonstration, the file name is "B4U_HUA_toy_genetics".

## Required software

The tools below will be needed for the proper execution of the analysis plan.
The tools have been distributed as Docker/Singularity images.

|Tool|Version|Webpage|Download|
|:----|:----|:----|:----|
|bcftools|1.20|[https://www.htslib.org](https://www.htslib.org)|[download](https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2)|
|htslib|1.20|[https://www.htslib.org](https://www.htslib.org)|[download](https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2)|
|VCFtools|0.1.16|[https://vcftools.github.io/](https://vcftools.github.io/)|[download](https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz)|
|PLINK|1.90|[https://www.cog-genomics.org/plink/](https://www.cog-genomics.org/plink/)|[download](https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip)|
|REGENIE|3.5|[https://rgcgithub.github.io/regenie/](https://rgcgithub.github.io/regenie/)|[download](https://github.com/rgcgithub/regenie/releases/download/v3.5/regenie_v3.5.gz_x86_64_Linux.zip)|
|GCTA|1.94.1|[https://yanglab.westlake.edu.cn/software/gcta/](https://yanglab.westlake.edu.cn/software/gcta/)|[download](https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.1-linux-kernel-3-x86_64.zip)|
|KING|2.3.2|[https://kingrelatedness.com/](https://kingrelatedness.com/)|[download](https://kingrelatedness.com/Linux-king.tar.gz)|
|MR-MEGA|0.2|[https://genomics.ut.ee/en/tools](https://genomics.ut.ee/en/tools)|[download](https://tools.gi.ut.ee/tools/MR-MEGA_v0.2.zip)|

# Proposed workflow

## Phenotype file data management

To install the packages from CRAN, please use

```r
Rscript \
-e '{
install.packages(c("tidyr", "dplyr", "lubridate", "data.table"))
}'
```

Two [source files](https://drive.google.com/drive/folders/1UPzrMvIjzdXUsTv9BI_DUbrBbhB6ffVi?usp=drive_link) need to be loaded in order to calculate (i) the zBMI score for samples aged < 18 years and (ii) the summary statistics for the phenotype of interest (BMI/weight). Then, the following steps are performed:

* Appropriate data transformation of the phenotype file to the required format.
* Generation of summary statistics for BMI and weight per age group (childhood/adolescence/elderly).
* Check for presence of heteroscedasticity of regression residuals. If no issues are detected, proceed by carrying out the modeling step with `fastGWA`.

```r
Rscript \
-e '{
###########################################################################
#
# Libraries
#
###########################################################################

lib <- c("tidyr",
         "dplyr",
         "lubridate",
         "data.table")
lapply(lib, require, character.only = TRUE)

###########################################################################
#
# Paths
#
###########################################################################

main_path <- "/path_to_folder/"

###########################################################################
#
# Source files
#
###########################################################################

source(paste0(main_path, 
              "zBMI_calculation", 
              ".R"))

source(paste0(main_path, 
              "gen_sumstats_functions", 
              ".R"))

###########################################################################
#
# Data management - phenotype file
#
###########################################################################

phen        <- data.table::fread(paste0(main_path, "B4U_HUA_toy_dataset.csv"), header = TRUE)
class(phen) <- class(phen)[2] 

#---- Keep only the fields of interest:
phen$FID_IID   <- paste(phen$FID, phen$IID, sep = "_")
phen_cols_keep <- c("FID_IID",
                    "FID",
                    "IID",
                    "date",
                    "sex",
                    "age",
                    "bw",
                    "bh",
                    "bmi")

phen <- phen[phen_cols_keep]

#---- Remove individuals with only a baseline record:
phen <-
  phen %>% 
  dplyr::group_by(FID_IID) %>% 
  filter(n() >= 2 )

#---- Deduce age at followup, when it is missing:
phen$date <- lubridate::dmy(phen$date)

# Assumes that sample-specific records are ordered by date: 
phen <-
  phen %>% 
  dplyr::group_by(FID_IID) %>% 
  dplyr::mutate(age_missing = round(interval(min(date), date) / years(1), 1 )) %>% 
  dplyr::mutate(age2 = ifelse( is.na(age), age[1] + cumsum(c(0, age_missing[-1])), age) ) %>% 
  dplyr::select(-c(age, age_missing)) %>% 
  dplyr::rename(age = age2)

#---- Calculate zBMI for kids/adolsescents:
zbmi      <- vector("numeric", length = nrow(phen))
for(i in 1:nrow(phen)) zbmi[i] <- calculate_zbmi(phen[i,'sex'], phen[i,'age'], phen[i,'bmi'])
phen$zbmi <- round( as.numeric(zbmi), 2)

#---- Subset by age a group:
phen_data_childhood_IDs <-
  phen %>% 
  dplyr::group_by(FID_IID) %>% 
  filter( age[1] <= 18) %>% 
  distinct(FID_IID) %>% 
  ungroup()

phen_data_adulthood_IDs <-
  phen %>% 
  dplyr::group_by(FID_IID) %>% 
  filter( age[1] > 18 & age[1] < 65) %>% 
  distinct(FID_IID) %>% 
  ungroup()

phen_data_elders_IDs <-
  phen %>% 
  dplyr::group_by(FID_IID) %>% 
  filter( age[1] >= 65) %>% 
  distinct(FID_IID) %>% 
  ungroup()

# Check:
length(unique(phen$FID_IID)) == (nrow(phen_data_childhood_IDs) + nrow(phen_data_adulthood_IDs) + nrow(phen_data_elders_IDs))

phen_data_childhood <- subset(phen, FID_IID %in% phen_data_childhood_IDs$FID_IID)
phen_data_adulthood <- subset(phen, FID_IID %in% phen_data_adulthood_IDs$FID_IID)
phen_data_elders    <- subset(phen, FID_IID %in% phen_data_elders_IDs$FID_IID)

# Check:
length(unique(phen$FID_IID)) == (
  length(unique(phen_data_childhood$FID_IID)) +
  length(unique(phen_data_adulthood$FID_IID)) + 
  length(unique(phen_data_elders$FID_IID)) )

#---- Convert to the desired data format:
phen_childhood_wide <- 
  phen_data_childhood %>%
  dplyr::group_by(FID_IID) %>%
  dplyr::summarise(
    Sex    = unique(sex),
    Age    = paste(age, collapse = ";"),
    Height = paste(bh,  collapse = ";"),
    Weight = paste(bw,  collapse = ";"),
    zBMI   = paste(zbmi, collapse = ";")
  )

phen_adulthood_wide <- 
  phen_data_adulthood %>%
  dplyr::group_by(FID_IID) %>%
  dplyr::summarise(
    Sex    = unique(sex),
    Age    = paste(age, collapse = ";"),
    Height = paste(bh,  collapse = ";"),
    Weight = paste(bw,  collapse = ";"),
    BMI    = paste(bmi, collapse = ";")
  )

phen_elders_wide <- 
  phen_data_elders %>%
  dplyr::group_by(FID_IID) %>%
  dplyr::summarise(
    Sex    = unique(sex),
    Age    = paste(age, collapse = ";"),
    Height = paste(bh,  collapse = ";"),
    Weight = paste(bw,  collapse = ";"),
    BMI    = paste(bmi, collapse = ";")
  )

#---- Clear up some memory:
rm("phen_cols_keep", 
   "phen_data_adulthood",    
   "phen_data_adulthood_IDs", 
   "phen_data_childhood",     
   "phen_data_childhood_IDs",
   "phen_data_elders",        
   "phen_data_elders_IDs", 
   "phen",
   "zbmi")

###########################################################################
#
# Generation of summary statistics
#
###########################################################################

#---- BMI:
sumstats_BMI_childhood <- 
   generate_sumstats(input_data = phen_childhood_wide,
                     pheno      = "BMI",       
                     age_group  = "Childhood") 

sumstats_BMI_adulthood <- 
  generate_sumstats(input_data = phen_adulthood_wide,
                    pheno      = "BMI",      
                    age_group  = "Adulthood") 

sumstats_BMI_elders <- 
  generate_sumstats(input_data = phen_elders_wide,
                    pheno      = "BMI",     
                    age_group  = "Elderly") 

#---- Weight:
sumstats_Weight_childhood <- 
  generate_sumstats(input_data = phen_childhood_wide,
                    pheno      = "Weight",   
                    age_group  = "Childhood") 

sumstats_Weight_adulthood <- 
  generate_sumstats(input_data = phen_adulthood_wide,
                    pheno      = "Weight",    
                    age_group  = "Adulthood") 

sumstats_Weight_elders <- 
  generate_sumstats(input_data = phen_elders_wide,
                    pheno      = "Weight",
                    age_group  = "Elderly")


# ΝΟΤΕ: In the HUA biobank, all samples who had a baseline measurement during childhood,
#       had a followup measurement during adulthood. Therefore, the proposed process does
#       not calculate the phenotype rate of change (slope).

write.table(sumstats_BMI_adulthood[sumstats_BMI_adulthood$N >= 2, c(1:2,7:ncol(sumstats_BMI_adulthood))],
            paste0(main_path, "Toy_bmi_summary_adults.txt"),
            sep      = "\t",
            quote    = F,
            row.names= F)

write.table(sumstats_BMI_elders[sumstats_BMI_elders$N >= 2, c(1:2,7:ncol(sumstats_BMI_elders))],
            paste0(main_path, "Toy_bmi_summary_elderly.txt"),
            sep      = "\t",
            quote    = F,
            row.names= F)

write.table(sumstats_Weight_adulthood[sumstats_Weight_adulthood$N >= 2, c(1:2,7:ncol(sumstats_Weight_adulthood))],
            paste0(main_path, "Toy_weight_summary_adults.txt"),
            sep      = "\t",
            quote    = F,
            row.names= F)

write.table(sumstats_Weight_elders[sumstats_Weight_elders$N >= 2, c(1:2,7:ncol(sumstats_Weight_elders))],
            paste0(main_path, "Toy_weight_summary_elderly.txt"),
            sep      = "\t",
            quote    = F,
            row.names= F)

#---- Export the unique sample IDs per age group:
write.table(sumstats_BMI_adulthood[sumstats_BMI_adulthood$N >= 2, 1],
            paste0(main_path, "Toy_bmi_adults_IDs.txt"),
            sep      = "\t",
            quote    = F,
            row.names= F,
            col.names= F)

write.table(sumstats_BMI_elders[sumstats_BMI_elders$N >= 2, 1],
            paste0(main_path, "Toy_bmi_elderly_IDs.txt"),
            sep      = "\t",
            quote    = F,
            row.names= F,
            col.names= F)

write.table(sumstats_Weight_adulthood[sumstats_Weight_adulthood$N >= 2, 1],
            paste0(main_path, "Toy_weight_adults_IDs.txt"),
            sep      = "\t",
            quote    = F,
            row.names= F,
            col.names= F)

write.table(sumstats_Weight_elders[sumstats_Weight_elders$N >= 2, 1],
            paste0(main_path, "Toy_weight_elderly_IDs.txt"),
            sep      = "\t",
            quote    = F,
            row.names= F,
            col.names= F)

###########################################################################
#
# Heteroscedasticity check
#
###########################################################################

check_sumstats_BMI_adulthood    <- sumstats_BMI_adulthood[sumstats_BMI_adulthood$N >= 2, ]
check_sumstats_BMI_elders       <- sumstats_BMI_elders[sumstats_BMI_elders$N >= 2, ]

check_sumstats_Weight_adulthood <- sumstats_Weight_adulthood[sumstats_Weight_adulthood$N >= 2, ]
check_sumstats_Weight_elders    <- sumstats_Weight_elders[sumstats_Weight_elders$N >= 2, ]

summary(check_sumstats_BMI_adulthood$S_Age)
summary(check_sumstats_BMI_elders$S_Age)

summary(check_sumstats_Weight_adulthood$S_Age)
summary(check_sumstats_Weight_elders$S_Age)
}'
```

**NOTE**

The toy example contains sample IDs whose baseline measurement was taken during childhood/adolescence and the followup measurement(s) were taken during adulthood. See object `phen_childhood_wide` above. During the process for the calculation of the summary statistics, these followup measurements would be discarded, leaving only one measurement per sample, not contributing to the rate of change of the desired phenotype. Therefore, these sample IDs are not analysed further in the context of the toy dataset.

## Principal Component projections

[Files](https://drive.google.com/drive/folders/1Y7jNcJmayhF9PdcmvBACFktSfFHcabWA) based on Principal Component Analysis of the current 1000 genomes data have been provided by the central analysis team. These files are:

- `pca_variants.txt`: a set of SNPs whose loadings are used to project the local genotypes to 1000 genomes and be used as covariates in subsequent GWAS.
- `loads_1000g.txt`: the aforementioned loadings.
- `means_1000g.txt`: means and standard deviations of the loadings required by
`flashpca`.

In this step the file `pca_variants.txt` are used in order to create a PLINK fileset with a subset of variants contained in this file. Then, these are used with the files `loads_1000g.txt` and `means_1000g.txt` along with the projection functionality of `flashpca` to create the PC covariates for your cohort.

**IMPORTANT**: Although the central analysis team will make sure that the  variants  used for PCA projections are common across all partners , it remains **your** responsibility to make sure that **all** the variants in the provided file `pca_variants.txt` are present in your cohort. If the cohort is not very small and has been imputed with HRC panel, this should not be a problem. Otherwise, please contact the central analysis team  **immediately**. It is also **your** responsibility to ensure that the major/minor alleles present in the PCA SNPs are in alignment with your dataset. A basic command is provided below to ensure this.

Firstly, extract the reference allele from `loads_1000g.txt`:

```bash
tail -n +2 loads_1000g.txt | cut -f 1,2 > refpos_1000g.txt
```

Then, create the PLINK files for PCA projection, switching some alleles if
necessary:

```bash
plink \
--bfile B4U_HUA_toy_genetics \
--a1-allele refpos_1000g.txt 2 1 \
--out B4U_HUA_toy_for_PCA \
--extract pca_variants.txt \
--make-bed
```

Then, project (use `flashpca_x86-64` instead of `flashpca` if not running fron
the Singularity image):

```
flashpca \
--bfile B4U_HUA_toy_for_PCA \
--inmeansd means_1000g.txt \
--inload loads_1000g.txt \
--project \
--outproj B4U_HUA_toy_projections.txt \
--verbose
```

## Regression modeling

### Regular time intervals - fastGWA

#### Preparation of covariate and phenotype files

```r
Rscript \
-e '{
###########################################################################
#
# Libraries
#
###########################################################################
lib <- c("tidyr",
		 "dplyr",
		 "lubridate",
		 "data.table")
lapply(lib, require, character.only = TRUE)

###########################################################################
#
# Paths
#
###########################################################################

main_path <- "//your_path//"
 
###########################################################################
#
# Data management - phenotype + PCA file
#
###########################################################################
proj           <- read.delim(paste0(main_path, "B4U_HUA_toy_projections.txt"))
proj$FID_IID   <- paste(proj$FID, proj$IID, sep = "_")
proj_cols_keep <- c("FID_IID",
                    "FID",
                    "IID",
                    paste0("PC",1:10))
proj           <- proj[proj_cols_keep]

# Round eigenvectors to the 3rd digit
for ( i in 4:ncol(proj )) proj[,i] <- round(proj[,i],3)

sumstats_BMI_adulthood <- 
 data.table::fread(paste0(main_path, "Toy_bmi_summary_adults.txt"),    
                   header = TRUE)
sumstats_BMI_elderly <- 
 data.table::fread(paste0(main_path, "Toy_bmi_summary_elderly.txt"),    
                   header = TRUE)
sumstats_weight_adulthood <- 
 data.table::fread(paste0(main_path, "Toy_weight_summary_adults.txt"),  
                   header = TRUE)
sumstats_weight_elderly   <- 
 data.table::fread(paste0(main_path, "Toy_weight_summary_elderly.txt"), 
                   header = TRUE)

class(sumstats_BMI_adulthood) <- 
 class(sumstats_BMI_elderly) <- 
  class(sumstats_weight_adulthood) <- 
   class(sumstats_weight_elderly) <- class(sumstats_weight_elderly)[2] 

sumstats_BMI_adulthood <- 
  sumstats_BMI_adulthood %>% 
  left_join(proj, c("FID_IID" = "FID_IID"))

sumstats_BMI_elderly <- 
  sumstats_BMI_elderly %>% 
  left_join(proj, c("FID_IID" = "FID_IID"))

sumstats_weight_adulthood <- 
  sumstats_weight_adulthood %>% 
  left_join(proj, c("FID_IID" = "FID_IID"))

sumstats_weight_elderly <- 
  sumstats_weight_elderly %>% 
  left_join(proj, c("FID_IID" = "FID_IID"))

#---- Export datasets for BMI:
bmi_keep_adulthood       <- c("FID", "IID", "BMI_Beta")
bmi_qcovs_keep_adulthood <- c("FID", "IID", "M_Age", paste0("PC",1:10))
bmi_ccovs_keep_adulthood <- c("FID", "IID", "Sex")

write.table(sumstats_BMI_adulthood[bmi_keep_adulthood],
			file      = paste0(main_path, "Toy_bmi_adults_gcta_phenotype.txt"),
			row.names = FALSE,
			quote     = FALSE)

write.table(sumstats_BMI_adulthood[, bmi_qcovs_keep_adulthood, drop = FALSE],
			file      = paste0(main_path, "Toy_bmi_adults_gcta_q_covariates.txt"),
			row.names = FALSE,
			quote     = FALSE)

write.table(sumstats_BMI_adulthood[,bmi_ccovs_keep_adulthood, drop = FALSE],
			file      = paste0(main_path, "Toy_bmi_adults_gcta_c_covariates.txt"),
			row.names = FALSE,
			quote     = FALSE)

write.table(sumstats_BMI_elderly[bmi_keep_adulthood],
			file      = paste0(main_path, "Toy_bmi_elderly_gcta_phenotype.txt"),
			row.names = FALSE,
			quote     = FALSE)

write.table(sumstats_BMI_elderly[, bmi_qcovs_keep_adulthood, drop = FALSE],
			file      = paste0(main_path, "Toy_bmi_elderly_gcta_q_covariates.txt"),
			row.names = FALSE,
			quote     = FALSE)

write.table(sumstats_BMI_elderly[,bmi_ccovs_keep_adulthood, drop = FALSE],
			file      = paste0(main_path, "Toy_bmi_elderly_gcta_c_covariates.txt"),
			row.names = FALSE,
			quote     = FALSE)

#---- Export datasets for weight:
weight_keep_adulthood       <- c("FID", "IID", "Weight_Beta")
weight_qcovs_keep_adulthood <- c("FID", "IID", "M_Age", paste0("PC",1:10))
weight_ccovs_keep_adulthood <- c("FID", "IID", "Sex")

write.table(sumstats_weight_adulthood[weight_keep_adulthood],
			file      = paste0(main_path, "Toy_weight_adults_gcta_phenotype.txt"),
			row.names = FALSE,
			quote     = FALSE)

write.table(sumstats_weight_adulthood[, weight_qcovs_keep_adulthood, drop = FALSE],
			file      = paste0(main_path, "Toy_weight_adults_gcta_q_covariates.txt"),
			row.names = FALSE,
			quote     = FALSE)

write.table(sumstats_weight_adulthood[,weight_ccovs_keep_adulthood, drop = FALSE],
			file      = paste0(main_path, "Toy_weight_adults_gcta_c_covariates.txt"),
			row.names = FALSE,
			quote     = FALSE)

write.table(sumstats_weight_elderly[weight_keep_adulthood],
			file      = paste0(main_path, "Toy_weight_elderly_gcta_phenotype.txt"),
			row.names = FALSE,
			quote     = FALSE)

write.table(sumstats_weight_elderly[, weight_qcovs_keep_adulthood, drop = FALSE],
			file      = paste0(main_path, "Toy_weight_elderly_gcta_q_covariates.txt"),
			row.names = FALSE,
			quote     = FALSE)

write.table(sumstats_weight_elderly[,weight_ccovs_keep_adulthood, drop = FALSE],
			file      = paste0(main_path, "Toy_weight_elderly_gcta_c_covariates.txt"),
			row.names = FALSE,
			quote     = FALSE)
	
#---- Export the unique sample IDs per age group:
write.table(sumstats_BMI_adulthood[c("FID", "IID")],
			paste0(main_path, "Toy_bmi_adults_IDs.txt"),
			sep      = "\t",
			quote    = F,
			row.names= F,
			col.names= F)

write.table(sumstats_BMI_elderly[c("FID", "IID")],
			paste0(main_path, "Toy_bmi_elderly_IDs.txt"),
			sep      = "\t",
			quote    = F,
			row.names= F,
			col.names= F)

write.table(sumstats_weight_adulthood[c("FID", "IID")],
			paste0(main_path, "Toy_weight_adults_IDs.txt"),
			sep      = "\t",
			quote    = F,
			row.names= F,
			col.names= F)

write.table(sumstats_weight_elderly[c("FID", "IID")],
			paste0(main_path, "Toy_weight_elderly_IDs.txt"),
			sep      = "\t",
			quote    = F,
			row.names= F,
			col.names= F)
}'
```

#### Preparation of the Genetic Relationships Matrix (GRM)

The GRM is a prerequisite for using GCTA and can be constructed with the GCTA software taking into account sample relatedness as calculated with the KING software. 

**Step 1: Find related individuals with KING**

Subset the `B4U_HUA_toy_genetics` dataset according to the age groups of interest.

```
plink \
--bfile B4U_HUA_toy_genetics \
--keep Toy_bmi_adults_IDs.txt \
--make-bed \
--out Toy_bmi_adults

plink \
--bfile B4U_HUA_toy_genetics \
--keep Toy_bmi_elderly_IDs.txt \
--make-bed \
--out Toy_bmi_elderly

plink \
--bfile B4U_HUA_toy_genetics \
--keep Toy_weight_adults_IDs.txt \
--make-bed \
--out Toy_weight_adults

plink \
--bfile B4U_HUA_toy_genetics \
--keep Toy_weight_elderly_IDs.txt \
--make-bed \
--out Toy_weight_elderly
```
**Step 2: Exclude monomorphic and low variance SNPs**

```
plink \
--bfile Toy_bmi_adults \
--out Toy_bmi_adults_f \
--make-bed \
--maf 0.01

plink \
--bfile Toy_bmi_elderly \
--out Toy_bmi_elderly_f \
--make-bed \
--maf 0.01

plink \
--bfile Toy_weight_adults \
--out Toy_weight_adults_f \
--make-bed \
--maf 0.01

plink \
--bfile Toy_weight_elderly \
--out Toy_weight_elderly_f \
--make-bed \
--maf 0.01
```

**Step 3: Find related individuals with KING**

```
king \
-b Toy_bmi_adults_f.bed \
--unrelated \
--degree 2 \
--prefix bmi_adults_individuals_

king \
-b Toy_bmi_elderly_f.bed \
--unrelated \
--degree 2 \
--prefix bmi_elderly_individuals_

king \
-b Toy_weight_adults_f.bed \
--unrelated \
--degree 2 \
--prefix weight_adults_individuals_

king \
-b Toy_weight_elderly_f.bed \
--unrelated \
--degree 2 \
--prefix weight_elderly_individuals_
```

**Step 4: Conduct LD-pruning excluding related samples**

```
plink \
--bfile Toy_bmi_adults_f \
--keep bmi_adults_individuals_unrelated.txt \
--indep 50 5 2 \
--out Toy_bmi_adults_pruned

plink \
--bfile Toy_bmi_elderly_f \
--keep bmi_elderly_individuals_unrelated.txt \
--indep 50 5 2 \
--out Toy_bmi_elderly_pruned

plink \
--bfile Toy_weight_adults_f \
--keep weight_adults_individuals_unrelated.txt \
--indep 50 5 2 \
--out Toy_weight_adults_pruned

plink \
--bfile Toy_weight_elderly_f \
--keep weight_elderly_individuals_unrelated.txt \
--indep 50 5 2 \
--out Toy_weight_elderly_pruned
```

**Step 5: Create PLINK files with the pruned variants and all samples**

```
plink \
--bfile Toy_bmi_adults_f \
--extract Toy_bmi_adults_pruned.prune.in \
--out Toy_bmi_adults_for_grm \
--make-bed

plink \
--bfile Toy_bmi_elderly_f \
--extract Toy_bmi_elderly_pruned.prune.in \
--out Toy_bmi_elderly_for_grm \
--make-bed

plink \
--bfile Toy_weight_adults_f \
--extract Toy_weight_adults_pruned.prune.in \
--out Toy_weight_adults_for_grm \
--make-bed

plink \
--bfile Toy_weight_elderly_f \
--extract Toy_weight_elderly_pruned.prune.in \
--out Toy_weight_elderly_for_grm \
--make-bed
```

**Step 6: Create GRM and GRM sparse matrices**

`--thread-num` can be set according to the system of each user. We set a small number for the toy dataset.

```
gcta64 \
--bfile Toy_bmi_adults_for_grm \
--make-grm \
--autosome \
--thread-num 4 \
--out Toy_bmi_adults_grm

gcta64 \
--grm Toy_bmi_adults_grm \
--make-bK-sparse 0.05 \
--out Toy_bmi_adults_grm_sparse

gcta64 \
--bfile Toy_bmi_elderly_for_grm \
--make-grm \
--autosome \
--thread-num 4 \
--out Toy_bmi_elderly_grm

gcta64 \
--grm Toy_bmi_elderly_grm \
--make-bK-sparse 0.05 \
--out Toy_bmi_elderly_grm_sparse

gcta64 \
--bfile Toy_weight_adults_for_grm \
--make-grm \
--autosome \
--thread-num 4 \
--out Toy_weight_adults_grm

gcta64 \
--grm Toy_weight_adults_grm \
--make-bK-sparse 0.05 \
--out Toy_weight_adults_grm_sparse

gcta64 \
--bfile Toy_weight_elderly_for_grm \
--make-grm \
--autosome \
--thread-num 4 \
--out Toy_weight_elderly_grm

gcta64 \
--grm Toy_weight_elderly_grm \
--make-bK-sparse 0.05 \
--out Toy_weight_elderly_grm_sparse
```

**Step 7: Perform pruning**

```
plink \
--bfile Toy_bmi_adults \
--out Toy_bmi_adults \
--indep-pairwise 50 5 0.2

plink \
--bfile Toy_bmi_elderly \
--out Toy_bmi_elderly \
--indep-pairwise 50 5 0.2

plink \
--bfile Toy_weight_adults \
--out Toy_weight_adults \
--indep-pairwise 50 5 0.2

plink \
--bfile Toy_weight_elderly \
--out Toy_weight_elderly \
--indep-pairwise 50 5 0.2
```

#### Execute GCTA fastGWA with the toy dataset

**Step 1: Keep sample IDs with complete covariate information**

```
plink \
--bfile Toy_bmi_adults \
--keep Toy_bmi_adults_IDs.txt \
--make-bed \
--out Toy_bmi_adults_2

plink \
--bfile Toy_bmi_elderly \
--keep Toy_bmi_elderly_IDs.txt \
--make-bed \
--out Toy_bmi_elderly_2

plink \
--bfile Toy_weight_adults \
--keep Toy_weight_adults_IDs.txt \
--make-bed \
--out Toy_weight_adults_2

plink \
--bfile Toy_weight_elderly \
--keep Toy_weight_elderly_IDs.txt \
--make-bed \
--out Toy_weight_elderly_2
```

**Step 2: Create the null model using pruned SNPs with GCTA**

```
plink \
--bfile Toy_bmi_adults_2 \
--out Toy_bmi_adults_null \
--exclude Toy_bmi_adults_pruned.prune.out \
--make-bed

plink \
--bfile Toy_bmi_elderly_2 \
--out Toy_bmi_elderly_null \
--exclude Toy_bmi_elderly_pruned.prune.out \
--make-bed

plink \
--bfile Toy_weight_adults_2 \
--out Toy_weight_adults_null \
--exclude Toy_weight_adults_pruned.prune.out \
--make-bed

plink \
--bfile Toy_weight_elderly_2 \
--out Toy_weight_elderly_null \
--exclude Toy_weight_elderly_pruned.prune.out \
--make-bed
```

**Step 3: Perform final GWAS by also using the null model with GCTA**

```
gcta64 \
--fastGWA-mlm \
--est-vg HE \
--bfile Toy_bmi_adults_null \
--grm-sparse Toy_bmi_adults_grm_sparse \
--pheno Toy_bmi_adults_gcta_phenotype.txt \
--qcovar Toy_bmi_adults_gcta_q_covariates.txt \
--covar Toy_bmi_adults_gcta_c_covariates.txt \
--model-only \
--thread-num 4 \
--out Toy_bmi_adults_fit_gcta

gcta64 \
--fastGWA-mlm \
--est-vg HE \
--bfile Toy_bmi_elderly_null \
--grm-sparse Toy_bmi_elderly_grm_sparse \
--pheno Toy_bmi_elderly_gcta_phenotype.txt \
--qcovar Toy_bmi_elderly_gcta_q_covariates.txt \
--covar Toy_bmi_elderly_gcta_c_covariates.txt \
--model-only \
--thread-num 4 \
--out Toy_bmi_elderly_fit_gcta

gcta64 \
--fastGWA-mlm \
--est-vg HE \
--bfile Toy_weight_adults_null \
--grm-sparse Toy_weight_adults_grm_sparse \
--pheno Toy_weight_adults_gcta_phenotype.txt \
--qcovar Toy_weight_adults_gcta_q_covariates.txt \
--covar Toy_weight_adults_gcta_c_covariates.txt \
--model-only \
--thread-num 4 \
--out Toy_weight_adults_fit_gcta

gcta64 \
--fastGWA-mlm \
--est-vg HE \
--bfile Toy_weight_elderly_null \
--grm-sparse Toy_weight_elderly_grm_sparse \
--pheno Toy_weight_elderly_gcta_phenotype.txt \
--qcovar Toy_weight_elderly_gcta_q_covariates.txt \
--covar Toy_weight_elderly_gcta_c_covariates.txt \
--model-only \
--thread-num 4 \
--out Toy_weight_elderly_fit_gcta
```

**Step 4: Perform final GWAS by also using the null model with GCTA** 

Cut-off low frequency variants with MAF = 0.01.

```
gcta64 \
--bfile Toy_bmi_adults_f \
--load-model Toy_bmi_adults_fit_gcta.fastGWA \
--maf 0.01 \
--out Toy_bmi_adults_fit_gcta_out \
--thread-num 4

gcta64 \
--bfile Toy_bmi_elderly_f \
--load-model Toy_bmi_elderly_fit_gcta.fastGWA \
--maf 0.01 \
--out Toy_bmi_elderly_fit_gcta_out \
--thread-num 4

gcta64 \
--bfile Toy_weight_adults_f \
--load-model Toy_weight_adults_fit_gcta.fastGWA \
--maf 0.01 \
--out Toy_weight_adults_fit_gcta_out \
--thread-num 4

gcta64 \
--bfile Toy_weight_elderly_f \
--load-model Toy_weight_elderly_fit_gcta.fastGWA \
--maf 0.01 \
--out Toy_weight_elderly_fit_gcta_out \
--thread-num 4
```

The summary statistics are in `.fastGWA`.

**Step 5: Compress summary statistics file and upload**

```
zip -r Toy_bmi_adults_fit_gcta_out.zip Toy_bmi_adults_fit_gcta_out.fastGWA
zip -r Toy_bmi_elderly_fit_gcta_out.zip Toy_bmi_elderly_fit_gcta_out.fastGWA

zip -r Toy_weight_adults_fit_gcta_out.zip Toy_weight_adults_fit_gcta_out.fastGWA
zip -r Toy_weight_elderly_fit_gcta_out.zip Toy_weight_elderly_fit_gcta_out.fastGWA
```

### Irregular time intervals - R

## Upload results

* Regular time intervals: please upload the zipped log files with the results of fitting the null model (`.fastGWA`) (containing estimates of heritability etc.).

* Irregular time intervals: The scripts will generate output files gwas.bmi.chr*.txt, gwas.bmi.chr*.err.txt, gwas.weight.chr*.txt and gwas.bmi.chr*.err.txt, where * are the numbers 1...22. Please zip together all these files and upload to the WP3 google drive in the **BETTER4U > WPS > WP > Weight_BMI_change > Results** folder. Please label the .zip file using the short name of your institution in the consortium listing. If using the R script for unrelated individuals, please report also the estimated ratio of the sampling noise to individual environmental noise (printed by the script).



