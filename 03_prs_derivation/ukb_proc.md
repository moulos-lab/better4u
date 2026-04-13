BETTER4U Tasks 3.7-3.8 - PRS evaluation with UKBB
================================================================================

## Authors

**Panagiotis Moulos** <sup>1,2</sup><br/>
<sup>1</sup>Harokopio University of Athens, <sup>2</sup>BSRC Alexander Fleming

### Contact details

* Panagiotis Moulos (moulos@fleming.gr)
* Jon Anders Eriksson (anders.eriksson@ut.ee)

# Introduction

This document outlines the validation of the several PRS versions created in the
frameworks of Tasks 3.7 and 3.8 of the BETTER4U WP3. It covers the analyses
performed in the case of PRS derivation with all participating cohorts as well
as evaluation within the LOO framework developed to provide genetic components
to the ML/AI processes within WP5.

# Evaluation process

## Available data

UKBB phenotypic data contain thousands of measurements across individuals,
visit centers and ancestries (mostly whie British). Detailed descriptions can be
found [here](https://biobank.ndph.ox.ac.uk/showcase).

The available data for our case are:

- UKB genotypic data in `BGEN+SAMPLE` Oxford formatPrs
- Phenotypic data in tabular format

After inspection of the phenotypic data, we end up using the following in our
validation for BMI:

- Biological sex
- Ethnic background at first visit
- Body mas index at first visit
- Age (and age<sup>2</sup>) at first visit
- The batch of genotype measurement
- 10 genetic principal components

In addition to the above, we use the following for phenotype filtering purposes:

- Genetic sex

Sex is coded as:

```bash
0 female
1 male
```

Ethnicity is coded as:

```bash
1 white
1001 British
1002 Irish
1003 Any other white
```

From the genotypic data, we are using imputed genetic data in `.bgen` format.

## Baseline BMI

### Phenotypic file preparation

We identified the following phenotype codings in the UKB `.tab` files:

| #column | tab id | phenotype |
| --- | --- | --- |
| 1 | f.eid | eid |
| 11 | f.31.0.0 | sex_f31_0_0 |
| 69 | f50.0.0 | standing_height_f50_0_0 |
| 70 | f50.1.0 | standing_height_f50_1_0 |
| 71 | f50.2.0 | standing_height_f50_2_0 |
| 72 | f50.3.0 | standing_height_f50_3_0 |
| 78 | f53.0.0 | date_of_attending_assessment_centre_f53_0_0 |
| 79 | f53.1.0 | date_of_attending_assessment_centre_f53_1_0 |
| 80 | f53.2.0 | date_of_attending_assessment_centre_f53_2_0 |
| 81 | f53.3.0 | date_of_attending_assessment_centre_f53_2_0 |
| 10957 | f.21000.0.0 | ethnic_background_f21000_0_0 |
| 10958 | f.21000.1.0 | ethnic_background_f21000_1_0 |
| 10959 | f.21000.2.0 | ethnic_background_f21000_2_0 |
| 10960 | f.21001.0.0 | body_mass_index_bmi_f21001_0_0 |
| 10961 | f.21001.1.0 | body_mass_index_bmi_f21001_1_0 |
| 10962 | f.21001.2.0 | body_mass_index_bmi_f21001_2_0 |
| 10963 | f.21001.3.0 | body_mass_index_bmi_f21001_3_0 |
| 10964 | f.21002.0.0 | weight_f21002_0_0 |
| 10965 | f.21002.1.0 | weight_f21002_1_0 |
| 10966 | f.21002.2.0 | weight_f21002_2_0 |
| 10967 | f.21002.3.0 | weight_f21002_3_0 |
| 10968 | f.21003.0.0 | age_when_attended_assessment_centre_f21003_0_0 |
| 11031 | f.22000.0.0 | genotype_measurement_batch_f22000_0_0 |
| 11032 | f.22001.0.0 | genetic_sex_f22001_0_0 |
| 11037 | f.22006.0.0 | genetic_ethnic_grouping_f22006_0_0 |
| 11040 | f.22009.0.1 | genetic_principal_components_f22009_0_1 |
| 11041 | f.22009.0.2 | genetic_principal_components_f22009_0_2 |
| 11042 | f.22009.0.3 | genetic_principal_components_f22009_0_3 |
| 11043 | f.22009.0.4 | genetic_principal_components_f22009_0_4 |
| 11044 | f.22009.0.5 | genetic_principal_components_f22009_0_5 |
| 11045 | f.22009.0.6 | genetic_principal_components_f22009_0_6 |
| 11046 | f.22009.0.7 | genetic_principal_components_f22009_0_7 |
| 11047 | f.22009.0.8 | genetic_principal_components_f22009_0_8 |
| 11048 | f.22009.0.9 | genetic_principal_components_f22009_0_9 |
| 11049 | f.22009.0.10 | genetic_principal_components_f22009_0_10 |
| 11050 | f.22009.0.11 | genetic_principal_components_f22009_0_11 |
| 11051 | f.22009.0.12 | genetic_principal_components_f22009_0_12 |
| 11052 | f.22009.0.13 | genetic_principal_components_f22009_0_13 |
| 11053 | f.22009.0.14 | genetic_principal_components_f22009_0_14 |
| 11054 | f.22009.0.15 | genetic_principal_components_f22009_0_15 |
| 11055 | f.22009.0.16 | genetic_principal_components_f22009_0_16 |
| 11056 | f.22009.0.17 | genetic_principal_components_f22009_0_17 |
| 11057 | f.22009.0.18 | genetic_principal_components_f22009_0_18 |
| 11058 | f.22009.0.19 | genetic_principal_components_f22009_0_19 |
| 11059 | f.22009.0.20 | genetic_principal_components_f22009_0_20 |

We can extract these measurements from the main `.tab` file:

```bash
cut -f1,11,69-72,78-81,10957-10967,10968,11031,11032,11037,11040-11059 ukb39726.tab > \
  ukb_fields_for_b4u_prs.txt
```

We then perform basic phenotypic QC in R:

```R
# UKBB phenos that we need
phe <- read.delim("ukb_fields_for_b4u_prs.txt.gz")

# Field name hash (based on the table we gathered)
field_hash <- c(
    "f.eid"="eid",
    "f.31.0.0"="sex_f31_0_0",
    "f.53.0.0"="date_of_attending_assessment_centre_f53_0_0",
    "f.53.1.0"="date_of_attending_assessment_centre_f53_1_0",
    "f.53.2.0"="date_of_attending_assessment_centre_f53_2_0",
    "f.53.3.0"="date_of_attending_assessment_centre_f53_3_0",
    "f.21000.0.0"="ethnic_background_f21000_0_0",
    "f.21000.1.0"="ethnic_background_f21000_1_0",
    "f.21000.2.0"="ethnic_background_f21000_2_0",
    "f.21001.0.0"="body_mass_index_bmi_f21001_0_0",
    "f.21001.1.0"="body_mass_index_bmi_f21001_1_0",
    "f.21001.2.0"="body_mass_index_bmi_f21001_2_0",
    "f.21001.3.0"="body_mass_index_bmi_f21001_3_0",
    "f.21002.0.0"="weight_f21002_0_0",
    "f.21002.1.0"="weight_f21002_1_0",
    "f.21002.2.0"="weight_f21002_2_0",
    "f.21002.3.0"="weight_f21002_3_0",
    "f.21003.0.0"="age_when_attended_assessment_centre_f21003_0_0",
    "f.22000.0.0"="genotype_measurement_batch_f22000_0_0",
    "f.22001.0.0"="genetic_sex_f22001_0_0",
    "f.22006.0.0"="genetic_ethnic_grouping_f22006_0_0",
    "f.22009.0.1"="genetic_principal_components_f22009_0_1",
    "f.22009.0.2"="genetic_principal_components_f22009_0_2",
    "f.22009.0.3"="genetic_principal_components_f22009_0_3",
    "f.22009.0.4"="genetic_principal_components_f22009_0_4",
    "f.22009.0.5"="genetic_principal_components_f22009_0_5",
    "f.22009.0.6"="genetic_principal_components_f22009_0_6",
    "f.22009.0.7"="genetic_principal_components_f22009_0_7",
    "f.22009.0.8"="genetic_principal_components_f22009_0_8",
    "f.22009.0.9"="genetic_principal_components_f22009_0_9",
    "f.22009.0.10"="genetic_principal_components_f22009_0_10",
    "f.22009.0.11"="genetic_principal_components_f22009_0_11",
    "f.22009.0.12"="genetic_principal_components_f22009_0_12",
    "f.22009.0.13"="genetic_principal_components_f22009_0_13",
    "f.22009.0.14"="genetic_principal_components_f22009_0_14",
    "f.22009.0.15"="genetic_principal_components_f22009_0_15",
    "f.22009.0.16"="genetic_principal_components_f22009_0_16",
    "f.22009.0.17"="genetic_principal_components_f22009_0_17",
    "f.22009.0.18"="genetic_principal_components_f22009_0_18",
    "f.22009.0.19"="genetic_principal_components_f22009_0_19",
    "f.22009.0.20"="genetic_principal_components_f22009_0_20"
)
names(phe) <- field_hash[names(phe)]

# The fields that we need for filtering
fields <- c(
    "eid",
    "sex_f31_0_0",
    "ethnic_background_f21000_0_0",
    "body_mass_index_bmi_f21001_0_0",
    "age_when_attended_assessment_centre_f21003_0_0",
    "genotype_measurement_batch_f22000_0_0",
    "genetic_sex_f22001_0_0",
    "genetic_principal_components_f22009_0_1",
    "genetic_principal_components_f22009_0_2",
    "genetic_principal_components_f22009_0_3",
    "genetic_principal_components_f22009_0_4",
    "genetic_principal_components_f22009_0_5",
    "genetic_principal_components_f22009_0_6",
    "genetic_principal_components_f22009_0_7",
    "genetic_principal_components_f22009_0_8",
    "genetic_principal_components_f22009_0_9",
    "genetic_principal_components_f22009_0_10",
    "genetic_principal_components_f22009_0_11",
    "genetic_principal_components_f22009_0_12",
    "genetic_principal_components_f22009_0_13",
    "genetic_principal_components_f22009_0_14",
    "genetic_principal_components_f22009_0_15",
    "genetic_principal_components_f22009_0_16",
    "genetic_principal_components_f22009_0_17",
    "genetic_principal_components_f22009_0_18",
    "genetic_principal_components_f22009_0_19",
    "genetic_principal_components_f22009_0_20"
)
phe <- phe[,fields,drop=FALSE]

# Remove missing data for this simple study
phe_complete <- complete.cases(phe)
phe <- phe[phe_complete,,drop=FALSE]

# Primary filters in phenos:
# Keep:
# phe$sex_f31_0_0 == phe$genetic_sex_f22001_0_0
# phe$ethnic_background_f21000_0_0==1001
filt_gen <- phe$ethnic_background_f21000_0_0==1001 & 
    phe$sex_f31_0_0 == phe$genetic_sex_f22001_0_0
# BMI outliers?
m_bmi <- median(phe$body_mass_index_bmi_f21001_0_0)
i_bmi <- IQR(phe$body_mass_index_bmi_f21001_0_0)
bmi_filt <- phe$body_mass_index_bmi_f21001_0_0 > m_bmi + 3*i_bmi |
    phe$body_mass_index_bmi_f21001_0_0 < m_bmi - 3*i_bmi
# Filter
filt <- filt_gen & !bmi_filt

# Final and also assign age^2
phe_keep <- phe[filt,-c(3,7)]
phe_keep$age_when_attended_assessment_centre_f21003_0_0_sq <- 
    phe_keep$age_when_attended_assessment_centre_f21003_0_0^2

write.table(phe_keep,file="covariates_20pc.txt",quote=FALSE,sep="\t",
    row.names=FALSE)

# Keep only first 10 instead of 20 PCs
phe_keep_10pc <- phe_keep[,-c(16:25)]
write.table(phe_keep_10pc,file="covariates_10pc.txt",quote=FALSE,sep="\t",
    row.names=FALSE)
```

### Assessment of Polygenic Risk Scores

We firstly define a workspace, where the UKB genotype files as well as the 
covariates file are available:

```bash
export WORKSPACE=/media/raid/collaborations/ukb
mkdir -p $WORKSPACE && cd $WORKSPACE

git clone https://github.com/moulos-lab/better4u.git
```

#### All-inclusive cohorts

As previously [described](https://github.com/moulos-lab/better4u/blob/main/03_prs_derivation/README.md#prs-processes), we have derived 5 versions of a PRS 
using all available participating cohorts. The files are:

- `b4u_bmi_prscs_original.prs`
- `b4u_bmi_prscs_robust.prs`
- `b4u_bmi_prscs_recalibrated.prs`
- `b4u_bmi_sbrc_tgp.prs`
- `b4u_bmi_sbrc_ukb.prs`

We can retrieve these files from Google Drive:

```bash
# b4u_bmi_prscs_original.prs
gdown 1CrCfcyP06Wv-BSG1wAsHJ-3fdkLI81yw --output b4u_bmi_prscs_original.prs

# b4u_bmi_prscs_robust.prs
gdown 1gbxXQzf4W_SM3iED6mE1nrjB9ffeUAdz --output b4u_bmi_prscs_robust.prs

# b4u_bmi_prscs_recalibrated.prs
gdown 1nXYLr4vY4oKq9QWFpaBa3gzb1Q-vyjD2 --output b4u_bmi_prscs_recalibrated.prs

# b4u_bmi_sbrc_tgp.prs
gdown 1BfEC6LiJevPnukgS6ozkPih_sVnbncJI --output b4u_bmi_sbrc_tgp.prs

# b4u_bmi_sbrc_ukb.prs
gdown 1PDn0mHBBo3nPKaK6iCvQuhjc5LKFK2qb --output b4u_bmi_sbrc_ukb.prs
``` 

Below we produce evaluation metrics for UKB after basic phenotypic QC. To 
achieve maximum PRS coverage we skip QC in genotypic data.

```R
WORKSPACE <- Sys.getenv("WORKSPACE")

source(file.path(WORKSPACE,"better4u","03_prs_derivation","evalfuns.R"))

trait="body_mass_index_bmi_f21001_0_0"
genoBase <- "ukb"
covFile <- "covariates.txt"

# PRS-CS with 1000 genomes LD - original version
prsFile_PRSCS_ORG <- "b4u_bmi_prscs_original.prs"
# We will convert to BIM just once
sanFile_PRSCS_ORG <- sanitizePrs(prsFile_PRSCS_ORG,genoBase,perChr=TRUE,
    bgen=TRUE,ukb=TRUE,permaBim=TRUE,from="ready")
M_PRSCS_ORG <- evalPrs(sanFile_PRSCS_ORG,covFile,trait,genoBase,iidCol=1,
    bgen=TRUE,perChr=TRUE,ukb=TRUE,remFile="removed.txt")

# PRS-CS with 1000 genomes LD - robust SNPs - original weights
prsFile_PRSCS_ROB <- "b4u_bmi_prscs_robust.prs"
# We have the BIMs, no need for bgen or ukb to be TRUE
sanFile_PRSCS_ROB <- sanitizePrs(prsFile_PRSCS_ROB,genoBase,perChr=TRUE,
    from="ready",rc=0.4)
M_PRSCS_ROB <- evalPrs(sanFile_PRSCS_ROB,covFile,trait,genoBase,iidCol=1,
    bgen=TRUE,perChr=TRUE,ukb=TRUE,remFile="removed.txt")

# PRS-CS with 1000 genomes LD - robust SNPs - original weights
prsFile_PRSCS_REC <- "b4u_bmi_prscs_recalibrated.prs"
sanFile_PRSCS_REC <- sanitizePrs(prsFile_PRSCS_REC,genoBase,perChr=TRUE,
    from="ready",rc=0.4)
M_PRSCS_REC <- evalPrs(sanFile_PRSCS_REC,covFile,trait,genoBase,iidCol=1,
    bgen=TRUE,perChr=TRUE,ukb=TRUE,remFile="removed.txt")

# SBayesRC with custom 1000 genomes LD
prsFile_SBRC_TGP <- "b4u_bmi_sbrc_tgp.prs"
sanFile_SBRC_TGP <- sanitizePrs(prsFile_SBRC_TGP,genoBase,perChr=TRUE,
    from="ready",rc=0.4)
M_SBRC_TGP <- evalPrs(sanFile_SBRC_TGP,covFile,trait,genoBase,iidCol=1,
    bgen=TRUE,perChr=TRUE,ukb=TRUE,remFile="removed.txt")

# SBayesRC with built-in UKB LD
prsFile_SBRC_UKB <- "b4u_bmi_sbrc_ukb.prs"
sanFile_SBRC_UKB <- sanitizePrs(prsFile_SBRC_UKB,genoBase,perChr=TRUE,
    from="ready",rc=0.4)
M_SBRC_UKB <- evalPrs(sanFile_SBRC_UKB,covFile,trait,genoBase,iidCol=1,
    bgen=TRUE,perChr=TRUE,ukb=TRUE,remFile="removed.txt")

# Some additional values regarding coverage
N_PRSCS_ORG <- as.numeric(countLines(prsFile_PRSCS_ORG)) - 1
N_PRSCS_ROB <- as.numeric(countLines(prsFile_PRSCS_ROB)) - 1
N_PRSCS_REC <- as.numeric(countLines(prsFile_PRSCS_REC)) - 1
N_SBRC_TGP <- as.numeric(countLines(prsFile_SBRC_TGP)) - 1
N_SBRC_UKB <- as.numeric(countLines(prsFile_SBRC_UKB)) - 1

# Gather metrics to share and PRS (just values, no ids of any kind)
metrics <- data.frame(
    original=formatC(M_PRSCS_ORG$metrics,digits=6),
    robust=formatC(M_PRSCS_ROB$metrics,digits=6),
    recalibrated=formatC(M_PRSCS_REC$metrics,digits=6),
    sbayes_tgp=formatC(M_SBRC_TGP$metrics,digits=6),
    sbayes_ukb=formatC(M_SBRC_UKB$metrics,digits=6)
)

# Additonal data
add <- as.data.frame(rbind(
    as.integer(c(N_PRSCS_ORG,N_PRSCS_ROB,N_PRSCS_REC,N_SBRC_TGP,N_SBRC_UKB)),
    round(100*as.numeric(metrics[nrow(metrics),])/
        c(N_PRSCS_ORG,N_PRSCS_ROB,N_PRSCS_REC,N_SBRC_TGP,N_SBRC_UKB),digits=2)
))
rownames(add) <- c("snps_total","coverage")
colnames(add) <- names(metrics)

# Final metrics data frame
finalMetrics <- rbind(metrics,add)

write.table(finalMetrics,file="b4u_bmi_prs_metrics_UKB.txt",sep="\t",
    quote=FALSE,col.names=NA)
```

#### LOO cohorts

We repeat the same process for score files produced after averaging betas for
the three methods applied in LOO which resulted in three PRS files with averaged
betas across 9 LOO runs. These files are:

- `b4u_bmi_prscs_tgp_loo_consensus.prs`
- `b4u_bmi_sbrc_tgp_loo_consensus.prs`
- `b4u_bmi_sbrc_ukb_loo_consensus.prs`

We can retrieve these files from Google Drive:

```bash
# b4u_bmi_prscs_tgp_loo_consensus.prs
gdown  --output b4u_bmi_prscs_tgp_loo_consensus.prs

# b4u_bmi_sbrc_tgp_loo_consensus.prs
gdown  --output b4u_bmi_sbrc_tgp_loo_consensus.prs

# b4u_bmi_sbrc_ukb_loo_consensus.prs
gdown  --output b4u_bmi_sbrc_ukb_loo_consensus.prs
```

Below we produce evaluation metrics for UKB after basic phenotypic QC. To 
achieve maximum PRS coverage we skip QC in genotypic data.

```R
WORKSPACE <- Sys.getenv("WORKSPACE")

source(file.path(WORKSPACE,"better4u","03_prs_derivation","evalfuns.R"))

trait="body_mass_index_bmi_f21001_0_0"
genoBase <- "ukb"
covFile <- "covariates.txt"
#covFile <- "covariates_10pc.txt"

# PRS-CS with 1000 genomes LD - original version
prsFile_PRSCS_TGP <- "b4u_bmi_prscs_tgp_loo_consensus.prs"
# We have the BIMs from the previous all-inclusive run
sanFile_PRSCS_TGP <- sanitizePrs(prsFile_PRSCS_TGP,genoBase,perChr=TRUE,
    from="ready",rc=0.4)
M_PRSCS_TGP <- evalPrs(sanFile_PRSCS_TGP,covFile,trait,genoBase,iidCol=1,
    bgen=TRUE,perChr=TRUE,ukb=TRUE,remFile="removed.txt")

# SBayesRC with custom 1000 genomes LD
prsFile_SBRC_TGP <- "b4u_bmi_sbrc_tgp_loo_consensus.prs"
sanFile_SBRC_TGP <- sanitizePrs(prsFile_SBRC_TGP,genoBase,perChr=TRUE,
    from="ready",rc=0.4)
M_SBRC_TGP <- evalPrs(sanFile_SBRC_TGP,covFile,trait,genoBase,iidCol=1,
    bgen=TRUE,perChr=TRUE,ukb=TRUE,remFile="removed.txt")

# SBayesRC with built-in UKB LD
prsFile_SBRC_UKB <- "b4u_bmi_sbrc_ukb_loo_consensus.prs"
sanFile_SBRC_UKB <- sanitizePrs(prsFile_SBRC_UKB,genoBase,perChr=TRUE,
    from="ready",rc=0.4)
M_SBRC_UKB <- evalPrs(sanFile_SBRC_UKB,covFile,trait,genoBase,iidCol=1,
    bgen=TRUE,perChr=TRUE,ukb=TRUE,remFile="removed.txt")

# Some additional values regarding coverage
N_PRSCS_TGP <- as.numeric(countLines(prsFile_PRSCS_ORG)) - 1
N_SBRC_TGP <- as.numeric(countLines(prsFile_SBRC_TGP)) - 1
N_SBRC_UKB <- as.numeric(countLines(prsFile_SBRC_UKB)) - 1

# Gather metrics to share and PRS (just values, no ids of any kind)
metrics <- data.frame(
    prscs_tgp=formatC(M_PRSCS_TGP$metrics,digits=6),
    sbayes_tgp=formatC(M_SBRC_TGP$metrics,digits=6),
    sbayes_ukb=formatC(M_SBRC_UKB$metrics,digits=6)
)

# Additonal data
add <- as.data.frame(rbind(
    as.integer(c(N_PRSCS_TGP,N_SBRC_TGP,N_SBRC_UKB)),
    round(100*as.numeric(metrics[nrow(metrics),])/
        c(N_PRSCS_TGP,N_SBRC_TGP,N_SBRC_UKB),digits=2)
))
rownames(add) <- c("snps_total","coverage")
colnames(add) <- names(metrics)

# Final metrics data frame
finalMetrics <- rbind(metrics,add)

write.table(finalMetrics,file="b4u_bmi_prs_metrics_UKB_LOO.txt",sep="\t",
    quote=FALSE,col.names=NA)
#write.table(finalMetrics,file="b4u_bmi_prs_metrics_UKB_LOO_10pc.txt",sep="\t",
#    quote=FALSE,col.names=NA)
```

## BMI change

The BMI change case is generally more complex as:

1) It can be extracted only from different visits in UKB
2) The samples are one order of magnitude less (~40k)
3) *de novo* QC and PCA calculation must take place in UKB genetic data

### Phenotypic file preparation

We are using the same base phenotypes file as with the baseline BMI case in
order to prepare and QC:

```R
# UKBB phenos that we need
phe <- read.delim("ukb_fields_for_b4u_prs.txt.gz")

# Field name hash (based on the table we gathered)
field_hash <- c(
    "f.eid"="eid",
    "f.31.0.0"="sex_f31_0_0",
    "f.53.0.0"="date_of_attending_assessment_centre_f53_0_0",
    "f.53.1.0"="date_of_attending_assessment_centre_f53_1_0",
    "f.53.2.0"="date_of_attending_assessment_centre_f53_2_0",
    "f.53.3.0"="date_of_attending_assessment_centre_f53_3_0",
    "f.21000.0.0"="ethnic_background_f21000_0_0",
    "f.21000.1.0"="ethnic_background_f21000_1_0",
    "f.21000.2.0"="ethnic_background_f21000_2_0",
    "f.21001.0.0"="body_mass_index_bmi_f21001_0_0",
    "f.21001.1.0"="body_mass_index_bmi_f21001_1_0",
    "f.21001.2.0"="body_mass_index_bmi_f21001_2_0",
    "f.21001.3.0"="body_mass_index_bmi_f21001_3_0",
    "f.21003.0.0"="age_when_attended_assessment_centre_f21003_0_0",
    "f.22000.0.0"="genotype_measurement_batch_f22000_0_0",
    "f.22001.0.0"="genetic_sex_f22001_0_0",
    "f.22006.0.0"="genetic_ethnic_grouping_f22006_0_0",
    "f.22009.0.1"="genetic_principal_components_f22009_0_1",
    "f.22009.0.2"="genetic_principal_components_f22009_0_2",
    "f.22009.0.3"="genetic_principal_components_f22009_0_3",
    "f.22009.0.4"="genetic_principal_components_f22009_0_4",
    "f.22009.0.5"="genetic_principal_components_f22009_0_5",
    "f.22009.0.6"="genetic_principal_components_f22009_0_6",
    "f.22009.0.7"="genetic_principal_components_f22009_0_7",
    "f.22009.0.8"="genetic_principal_components_f22009_0_8",
    "f.22009.0.9"="genetic_principal_components_f22009_0_9",
    "f.22009.0.10"="genetic_principal_components_f22009_0_10",
    "f.22009.0.11"="genetic_principal_components_f22009_0_11",
    "f.22009.0.12"="genetic_principal_components_f22009_0_12",
    "f.22009.0.13"="genetic_principal_components_f22009_0_13",
    "f.22009.0.14"="genetic_principal_components_f22009_0_14",
    "f.22009.0.15"="genetic_principal_components_f22009_0_15",
    "f.22009.0.16"="genetic_principal_components_f22009_0_16",
    "f.22009.0.17"="genetic_principal_components_f22009_0_17",
    "f.22009.0.18"="genetic_principal_components_f22009_0_18",
    "f.22009.0.19"="genetic_principal_components_f22009_0_19",
    "f.22009.0.20"="genetic_principal_components_f22009_0_20"
)
names(phe) <- field_hash[names(phe)]

# The fields that we need for preprocessing (no PCs this time)
fields <- c(
    "eid",
    "sex_f31_0_0",
    "date_of_attending_assessment_centre_f53_0_0",
    "date_of_attending_assessment_centre_f53_1_0",
    "date_of_attending_assessment_centre_f53_2_0",
    "date_of_attending_assessment_centre_f53_3_0",
    "ethnic_background_f21000_0_0",
    "body_mass_index_bmi_f21001_0_0",
    "body_mass_index_bmi_f21001_1_0",
    "body_mass_index_bmi_f21001_2_0",
    "body_mass_index_bmi_f21001_3_0",
    "age_when_attended_assessment_centre_f21003_0_0",
    "genotype_measurement_batch_f22000_0_0",
    "genetic_sex_f22001_0_0"
)
phe <- phe[,fields,drop=FALSE]

# Define the number of samples for BMI change case
dbmi <- phe[,c(1,8:11)]
# Which visit has the most recorded cases?
nas <- apply(dbmi,2,function(x) length(which(!is.na(x))))
#                           eid body_mass_index_bmi_f21001_0_0
#                        502510                         499405
#body_mass_index_bmi_f21001_1_0 body_mass_index_bmi_f21001_2_0
#                         20299                          41692
#body_mass_index_bmi_f21001_3_0
#                           665
## Is seems that visit 2 is the most prominent
dbmi13 <- dbmi[,c(1,2,4)]
dbmi13c <- dbmi13[complete.cases(dbmi13),]
## 41625 cases, not bad!

# Define the data to be further filtered
rownames(dbmi13c) <- dbmi13c$eid
rownames(phe) <- phe$eid
phe <- phe[rownames(dbmi13c),,drop=FALSE]
# Still some sex missing
tmp <- phe[,c(2,7,12,13,14)]
cstmp <- complete.cases(tmp)
phe <- phe[cstmp,]

# Primary filters in phenos:
# Keep:
# phe$sex_f31_0_0 == phe$genetic_sex_f22001_0_0
# phe$ethnic_background_f21000_0_0==1001
filt_gen <- phe$ethnic_background_f21000_0_0==1001 & 
    phe$sex_f31_0_0 == phe$genetic_sex_f22001_0_0
# dBMI outliers?
dibmi <- phe$body_mass_index_bmi_f21001_2_0 - phe$body_mass_index_bmi_f21001_0_0
m_dibmi <- median(dibmi)
i_dibmi <- IQR(dibmi)
dibmi_filt <- dibmi > m_dibmi + 3*i_dibmi | dibmi < m_dibmi - 3*i_dibmi
# Filter
filt <- filt_gen & !dibmi_filt

# Final and also assign age^2, drop PCs, we won't need them
phe$delta_bmi = dibmi
#phe$age_when_attended_assessment_centre_f21003_0_0_sq <- 
#    phe$age_when_attended_assessment_centre_f21003_0_0^2
# For delta BMI
phe_keep <- phe[filt,-c(3:11,14)]
write.table(phe_keep,file="covariates_dbmi.txt",quote=FALSE,sep="\t",
    row.names=FALSE)
# For dBMI betas
phe_keep_beta <- phe[filt,-c(4,6,7,9,11,14:15)]
write.table(phe_keep_beta,file="covariates_beta_dbmi.txt",quote=FALSE,sep="\t",
    row.names=FALSE)
```

### Preparation of the regression response for δΒΜΙ

As the reponse variable when costruction the δBMI PRS from summary statistics
was consisted of the BMI change slopes (betas) instead of actual delta values,
we below perform the slope estimation for the UKB phenotypes prepared above
following the process and scripts distributed by A. Errikson:

```R
library(tidyverse)
library(lubridate)

d  <- read.delim("covariates_beta_dbmi.txt")
n <- nrow(d)

d$N <- d$M_Age <- d$S_Age <- d$M_BMI <- d$S_BMI <- d$BMI_Alpha <- 
    d$BMI_Beta <- rep(0,)0*(1:n);
 
d <- d %>%
    mutate(
        BMI=paste(body_mass_index_bmi_f21001_0_0,
            body_mass_index_bmi_f21001_2_0,sep=";")
    )

d <- d %>%
    mutate(
        # date_reading_r
        date_0 = ymd(date_of_attending_assessment_centre_f53_0_0),
        date_2 = ymd(date_of_attending_assessment_centre_f53_2_0),

        # difference in years 
        # date_2-date_0 to get days, then divide by 365.25
        years_diff = as.numeric(difftime(date_2,date_0,units="days"))/365.25,

        # Add the time difference to the original age
        age_when_attended_assessment_centre_f21003_2_0 = round(
            age_when_attended_assessment_centre_f21003_0_0 + years_diff,1
        ),
        
        # BMI1;BMI2
        BMI = paste(body_mass_index_bmi_f21001_0_0,
            body_mass_index_bmi_f21001_2_0,sep = ";"),

        # Age1;Age2
        Age = paste(age_when_attended_assessment_centre_f21003_0_0,
            age_when_attended_assessment_centre_f21003_2_0,sep = ";")
    )

for (k in 1:n) {
    # split age and BMI strings into vectors
    t <- as.numeric(strsplit(d$Age[k],";")[[1]])
    bmi <- suppressWarnings(as.numeric(strsplit(d$BMI[k],";")[[1]]))

    # select time points with BMI information and age in the right interval
    use <- !is.na(t) & !is.na(bmi) & t > 18

    # sub-set data to eligible time points
    t <- t[use] - 18 # offset so that Alpha is estimate of BMI at age 18
    bmi <- bmi[use]

    d$N[k] <- length(t)
    if (d$N[k] >= 2) {
        d$M_Age[k] <- mean(t);
        d$M_BMI[k] <- mean(bmi);
        # unbiased estimator of standard error of age in the sample
        d$S_Age[k] <- sqrt(sum((t-mean(t))^2)/(d$N[k]-1));
        # unbiased estimator of standard error of BMI in the sample
        d$S_BMI[k] <- sqrt(sum((bmi-mean(bmi))^2)/(d$N[k]-1));
        f <- summary(lm(bmi~t))
        # estimator of intercept at age 18
        d$BMI_Alpha[k] <- f$coefficients[1,1]
        # estimator of slope of weight
        d$BMI_Beta[k] <- f$coefficients[2,1]
    }
}

# Save data to file for individuals with >= 2 measurement points
write.table(d[d$N>=2,],"covariates_beta_dbmi_with_slopes.txt",sep="\t",
    quote=FALSE,row.names=F)
```

### Preparation of UKB data with δBMI samples

#### Conversion to PLINK 1.90 format

First, get a list of sample to be kept:

```bash
awk 'NR>1 {print $1" "$1}' covariates_beta_dbmi_with_slopes.txt > dbmi_keep.txt
```

and then perform the conversion while keeping the samples with BMI difference
available at the same time. `--memory` and `--threads` are adjusted for a 
machine with 64 available threads and 512GB RAM so as to process 22 chromosomes
concurrently. Adjust according to your environment.

```bash
for CHR in `seq 1 22`
do
  plink2 \
    --bgen ukb_chr${CHR}.bgen ref-first \
    --sample ukb_chr${CHR}.sample \
    --keep dbmi_keep.txt \
    --snps-only \
    --max-alleles 2 \
    --rm-dup force-first \
    --make-bed \
    --memory 22000 \
    --threads 3 \
    --out ukb_dbmi_chr${CHR} &
done
```

```bash
for CHR in `seq 1 22`
do
  echo ukb_dbmi_chr${CHR} >> merge.list
done

plink \
  --merge-list merge.list \
  --make-bed \
  --out ukb_dbmi
```

#### Basic SNP and sample filtering

Perform filtering with PLINK 1.9 to remove rare variants as well as lower sample 
and call rates. The following variant filters are recommended:

1. Variant call rate: >98%
2. Minor Allele Frequency: >0.05
3. Hardy-Weinberg equilibrium p-value: >10<sup>-6</sup>

The following sample filters are recommended:

1. Sample call rate: > 95%
2. Heterozygosity: median(heterozygosity) &plusmn; 3 &times; IQR
3. Identity By Descent: > 0.5 (optional)

Below, we sequentially apply variant and sample filters according to widely
accepted [best practices](https://onlinelibrary.wiley.com/doi/10.1002/sim.6605).

First, we calculate heterozygosities:

```bash
plink \
  --bfile ukb_dbmi \
  --out ukb_dbmi --het
```

Then, we create a file with sample names with heterozygosities within the limits
of our filter (sample filter #2 above):

```r
Rscript \
  -e '{
    hetData <- read.table("ukb_dbmi.het",header=TRUE,
        check.names=FALSE)
    rownames(hetData) <- hetData$IID
    heterozygosity <- 1 - hetData[,3]/hetData[,5]
    names(heterozygosity) <- rownames(hetData)
    avg <- median(heterozygosity,na.rm=TRUE)
    dev <- IQR(heterozygosity,na.rm=TRUE)
    keep <- heterozygosity > avg - 3*dev & heterozygosity < avg + 3*dev
    goodHet <- rownames(hetData)[keep]
    write.table(hetData[keep,c("FID","IID"),drop=FALSE],
      file="het_samples_pass.txt",col.names=FALSE,row.names=FALSE,quote=FALSE)
  }'
```

The file `het_samples_pass.txt` will be used after the next filters to compile 
a final list of samples to be kept in later analysis.

The following `plink` command will apply variant filters #1,2,3 and sample 
filter #1:

```bash
plink \
  --bfile ukb_dbmi \
  --out ukb_dbmi_tmp \
  --make-bed \
  --geno 0.02 \
  --maf 0.05 \
  --hwe 0.000001 \
  --mind 0.05
```

We now create files with variants and samples to *keep* (samples to keep are 
merged with those passing heterozygosity filters):

```r
cut -d" " -f1-2 ukb_dbmi_tmp.fam > generic_samples_pass.txt
cut -f2 ukb_dbmi_tmp.bim > generic_variants_pass.txt

Rscript \
  -e '{
    het <- read.table("het_samples_pass.txt")
    rownames(het) <- het[,2]
    gen <- read.table("generic_samples_pass.txt")
    rownames(gen) <- gen[,2]
    pass <- intersect(rownames(het),rownames(gen))
    write.table(gen[pass,,drop=FALSE],file="all_samples_pass.txt",
      col.names=FALSE,row.names=FALSE,quote=FALSE)
  }'
```

Based on the variant and sample content of the files `generic_variants_pass.txt`
and `all_samples_pass.txt` we create filtered PLINK files. These will be used
for kinship analysis filtering (can be skipped if not necesseary) and PCA.

```bash
plink \
  --bfile ukb_dbmi \
  --out ukb_dbmi_filtered \
  --extract generic_variants_pass.txt \
  --keep all_samples_pass.txt \
  --make-bed
```

Kinship analysis is performed with KING (no LD pruning):

```bash
king \
  -b ukb_dbmi_filtered.bed \
  --unrelated \
  --degree 2 \
  --prefix individuals_
```

The aforementioned command will create the file `individuals_unrelated.txt` 
which we use with PLINK to create the dataset with related individuals removed
(optional).

```bash
plink \
  --bfile ukb_dbmi_filtered \
  --keep individuals_unrelated.txt \
  --out ukb_dbmi_final \
  --make-bed
```

#### Perform LD pruning only δBMI data

Prior to PCA, it is common practice to not use variants in LD, so LD pruning is
performed.

* `--indep-pairwise [window size] [step size/variant count)] [R2]` controls the 
pruning parameters. e.g. indep `50 5 0.9` and generates a list of markers in 
approximate linkage equilibrium - takes 50 SNPs at a time and then shifts by 5 
for the window. R2 is the cut-off for linkage disequilibrium.
* We follow an iterative procedure similar to the
[FinnGen](https://finngen.gitbook.io/documentation/methods/phewas/quality-checks)
consortium, where R2 is decreased by a constant step until ~200,000 variants are
included in the pruned set.

```bash
mkdir pruned

TOTAL=10000000
R2=0.95
STEP=0.05

while [ $TOTAL -gt 200000 ]
do
  R2=$(echo "$R2-$STEP" | bc)
  
  plink \
    --bfile ukb_dbmi_final \
    --indep-pairwise 500 50 $R2 \
    --out ./pruned/ukb_dbmi_final
    
  TOTAL=$(wc -l ./pruned/*.in | awk '{print $1}')
  
  echo "Parameters --indep-pairwise 500 50 $R2 yield $TOTAL SNPs" >> \
    pruning.log
done
```

The final value of `$R2` is `0.3` and the number of SNPs is 180,309, so we
perform the final pruning:

```bash
rm -r ./pruned/*

FINAL_R2=`tail -1 pruning.log | awk '{print $(NF-3)}' | sed 's/^\./0./'`
plink \
  --bfile ukb_dbmi_final \
  --indep-pairwise 500 50 $FINAL_R2 \
  --out ukb_dbmi_final_pruned

plink \
  --bfile ukb_dbmi_final \
  --extract ukb_dbmi_final_pruned.prune.in \
  --make-bed \
  --out ukb_dbmi_final_pruned

rm *.nosex *.log
```

#### Principal Component Analysis with FlashPca2

Perform PCA using the FlashPCA2:

```
flashpca \
  --bfile ukb_dbmi_final_pruned \
  --ndim 20 \
  --outpc pcs_ukb_dbmi.txt \
  --outvec eig_ukb_dbmi.txt \
  --outload loads_ukb_dbmi.txt
```

#### Merge PCs with covariates file

The process is completed in R:

```R
covars <- read.delim("covariates_beta_dbmi_with_slopes.txt")
rownames(covars) <- covars$eid
covars <- covars[,-c(3:7,9:12,14:20)]
pcs <- read.delim("pcs_ukb_dbmi.txt")
rownames(pcs) <- pcs$IID
covars <- covars[rownames(pcs),,drop=FALSE]
covars20 <- cbind(covars,pcs[,3:ncol(pcs)])
covars10 <- cbind(covars,pcs[,3:12])
write.table(covars20,file="covariates_dbmi_20pc.txt",sep="\t",quote=FALSE,
    row.names=FALSE)
write.table(covars10,file="covariates_dbmi_10pc.txt",sep="\t",quote=FALSE,
    row.names=FALSE)
```

### Assessment of Polygenic Risk Scores for δΒΜΙ

We assume the same workspace as before.

#### All-inclusive cohorts

As previously [described](https://github.com/moulos-lab/better4u/tree/main/03_prs_derivation#weight-change), we have derived 5 versions of a PRS using all available participating cohorts. The 
files are:

- `b4u_wc_prscs_original.prs`
- `b4u_wc_prscs_robust.prs`
- `b4u_wc_prscs_recalibrated.prs`
- `b4u_wc_sbrc_tgp.prs`
- `b4u_wc_sbrc_ukb.prs`

We can retrieve these files from Google Drive:

```bash
# b4u_wc_prscs_original.prs
gdown 15136PvuqV9hpKV-Bc5LZTGv_FUv_Ew1H --output b4u_wc_prscs_original.prs

# b4u_wc_prscs_robust.prs
gdown 1c2OmdJJKtgvG8-gzwSHZLUHc0A0ad-hW --output b4u_wc_prscs_robust.prs

# b4u_wc_prscs_recalibrated.prs
gdown 17rrvwqluPgnl_fQPeIgFNsuwkFs9SaHB --output b4u_wc_prscs_recalibrated.prs

# b4u_wc_sbrc_tgp.prs
gdown 1LY1OZOdSCM3eyfMGxeQpK2JW1j40gAv1 --output b4u_wc_sbrc_tgp.prs

# b4u_wc_sbrc_ukb.prs
gdown 1hFGxhncwmY4lTMAOtS7wa5Hx7cOyngsC --output b4u_wc_sbrc_ukb.prs
``` 

Below we produce evaluation metrics for UKB after basic phenotypic QC. To 
achieve maximum PRS coverage we skip QC in genotypic data.

```R
WORKSPACE <- Sys.getenv("WORKSPACE")

source(file.path(WORKSPACE,"better4u","03_prs_derivation","evalfuns.R"))

trait <- "BMI_Beta"
genoBase <- "ukb_dbmi"
covFile <- "covariates_dbmi_20pcs_onlybeta.txt"
#covFile <- "covariates_dbmi_10pcs.txt"

# PRS-CS with 1000 genomes LD - original version
prsFile_PRSCS_ORG <- "b4u_wc_prscs_original.prs"
# We already have BIMs from the preprocessing
#sanFile_PRSCS_ORG <- sanitizePrs(prsFile_PRSCS_ORG,genoBase,perChr=TRUE,
#    from="ready",rc=0.4)
sanFile_PRSCS_ORG <- "b4u_wc_prscs_original.prs.san"
M_PRSCS_ORG <- evalPrs(sanFile_PRSCS_ORG,covFile,trait,genoBase,iidCol=1,
    perChr=TRUE,remFile="removed.txt",rc=0.4)

# PRS-CS with 1000 genomes LD - robust SNPs - original weights
prsFile_PRSCS_ROB <- "b4u_wc_prscs_robust.prs"
# We have the BIMs, no need for bgen or ukb to be TRUE
#sanFile_PRSCS_ROB <- sanitizePrs(prsFile_PRSCS_ROB,genoBase,perChr=TRUE,
#    from="ready",rc=0.4)
sanFile_PRSCS_ROB <- "b4u_wc_prscs_robust.prs.san"
M_PRSCS_ROB <- evalPrs(sanFile_PRSCS_ROB,covFile,trait,genoBase,iidCol=1,
    perChr=TRUE,remFile="removed.txt",rc=0.4)

# PRS-CS with 1000 genomes LD - robust SNPs - original weights
prsFile_PRSCS_REC <- "b4u_wc_prscs_recalibrated.prs"
#sanFile_PRSCS_REC <- sanitizePrs(prsFile_PRSCS_REC,genoBase,perChr=TRUE,
#    from="ready",rc=0.4)
sanFile_PRSCS_REC <- "b4u_wc_prscs_recalibrated.prs"
M_PRSCS_REC <- evalPrs(sanFile_PRSCS_REC,covFile,trait,genoBase,iidCol=1,
    perChr=TRUE,remFile="removed.txt",rc=0.4)

# SBayesRC with custom 1000 genomes LD
prsFile_SBRC_TGP <- "b4u_wc_sbrc_tgp.prs"
#sanFile_SBRC_TGP <- sanitizePrs(prsFile_SBRC_TGP,genoBase,perChr=TRUE,
#    from="ready",rc=0.4)
sanFile_SBRC_TGP <- "b4u_wc_sbrc_tgp.prs.san"
M_SBRC_TGP <- evalPrs(sanFile_SBRC_TGP,covFile,trait,genoBase,iidCol=1,
    perChr=TRUE,remFile="removed.txt",rc=0.4)

# SBayesRC with built-in UKB LD
prsFile_SBRC_UKB <- "b4u_wc_sbrc_ukb.prs"
#sanFile_SBRC_UKB <- sanitizePrs(prsFile_SBRC_UKB,genoBase,perChr=TRUE,
#    from="ready",rc=0.4)
sanFile_SBRC_UKB <- "b4u_wc_sbrc_ukb.prs.san"
M_SBRC_UKB <- evalPrs(sanFile_SBRC_UKB,covFile,trait,genoBase,iidCol=1,
    perChr=TRUE,remFile="removed.txt",rc=0.4)
    
save(M_PRSCS_ORG,M_PRSCS_ROB,M_PRSCS_REC,M_SBRC_TGP,M_SBRC_UKB,
    file="metrics_nopc_onlybeta.rda")
#save(M_PRSCS_ORG,M_PRSCS_ROB,M_PRSCS_REC,M_SBRC_TGP,M_SBRC_UKB,
#    file="metrics_10pc.rda")

# Some additional values regarding coverage
N_PRSCS_ORG <- as.numeric(countLines(prsFile_PRSCS_ORG)) - 1
N_PRSCS_ROB <- as.numeric(countLines(prsFile_PRSCS_ROB)) - 1
N_PRSCS_REC <- as.numeric(countLines(prsFile_PRSCS_REC)) - 1
N_SBRC_TGP <- as.numeric(countLines(prsFile_SBRC_TGP)) - 1
N_SBRC_UKB <- as.numeric(countLines(prsFile_SBRC_UKB)) - 1

# Gather metrics to share and PRS (just values, no ids of any kind)
metrics <- data.frame(
    original=formatC(M_PRSCS_ORG$metrics,digits=6),
    robust=formatC(M_PRSCS_ROB$metrics,digits=6),
    recalibrated=formatC(M_PRSCS_REC$metrics,digits=6),
    sbayes_tgp=formatC(M_SBRC_TGP$metrics,digits=6),
    sbayes_ukb=formatC(M_SBRC_UKB$metrics,digits=6)
)

# Additonal data
add <- as.data.frame(rbind(
    as.integer(c(N_PRSCS_ORG,N_PRSCS_ROB,N_PRSCS_REC,N_SBRC_TGP,N_SBRC_UKB)),
    round(100*as.numeric(metrics[nrow(metrics),])/
        c(N_PRSCS_ORG,N_PRSCS_ROB,N_PRSCS_REC,N_SBRC_TGP,N_SBRC_UKB),digits=2)
))
rownames(add) <- c("snps_total","coverage")
colnames(add) <- names(metrics)

# Final metrics data frame
finalMetrics <- rbind(metrics,add)

write.table(finalMetrics,file="b4u_wc_prs_metrics_UKB_nopc_onlybeta.txt",sep="\t",
    quote=FALSE,col.names=NA)
#write.table(finalMetrics,file="b4u_wc_prs_metrics_UKB_10pc.txt",sep="\t",
#    quote=FALSE,col.names=NA)
```

#### LOO cohorts

*WIP*

## Weight change

As with BMI change, the weight change case is generally more complex as:

1) It can be extracted only from different visits in UKB
2) The samples are one order of magnitude less (~40k)
3) *de novo* QC and PCA calculation must take place in UKB genetic data

### Phenotypic file preparation

We are using the same base phenotypes file as with the baseline BMI case in
order to prepare and QC:

```R
library(dplyr)
library(lubridate)

# UKBB phenos that we need
phe <- read.delim("ukb_fields_for_b4u_prs.txt.gz")

# Field name hash (based on the table we gathered)
field_hash <- c(
    "f.eid"="eid",
    "f.31.0.0"="sex_f31_0_0",
    "f.50.0.0"="standing_height_f50_0_0",
    "f.50.1.0"="standing_height_f50_1_0",
    "f.50.2.0"="standing_height_f50_2_0",
    "f.50.3.0"="standing_height_f50_3_0",
    "f.53.0.0"="date_of_attending_assessment_centre_f53_0_0",
    "f.53.1.0"="date_of_attending_assessment_centre_f53_1_0",
    "f.53.2.0"="date_of_attending_assessment_centre_f53_2_0",
    "f.53.3.0"="date_of_attending_assessment_centre_f53_3_0",
    "f.21000.0.0"="ethnic_background_f21000_0_0",
    "f.21000.1.0"="ethnic_background_f21000_1_0",
    "f.21000.2.0"="ethnic_background_f21000_2_0",
    "f.21001.0.0"="body_mass_index_bmi_f21001_0_0",
    "f.21001.1.0"="body_mass_index_bmi_f21001_1_0",
    "f.21001.2.0"="body_mass_index_bmi_f21001_2_0",
    "f.21001.3.0"="body_mass_index_bmi_f21001_3_0",
    "f.21002.0.0"="weight_f21002_0_0",
    "f.21002.1.0"="weight_f21002_1_0",
    "f.21002.2.0"="weight_f21002_2_0",
    "f.21002.3.0"="weight_f21002_3_0",
    "f.21003.0.0"="age_when_attended_assessment_centre_f21003_0_0",
    "f.22000.0.0"="genotype_measurement_batch_f22000_0_0",
    "f.22001.0.0"="genetic_sex_f22001_0_0",
    "f.22006.0.0"="genetic_ethnic_grouping_f22006_0_0",
    "f.22009.0.1"="genetic_principal_components_f22009_0_1",
    "f.22009.0.2"="genetic_principal_components_f22009_0_2",
    "f.22009.0.3"="genetic_principal_components_f22009_0_3",
    "f.22009.0.4"="genetic_principal_components_f22009_0_4",
    "f.22009.0.5"="genetic_principal_components_f22009_0_5",
    "f.22009.0.6"="genetic_principal_components_f22009_0_6",
    "f.22009.0.7"="genetic_principal_components_f22009_0_7",
    "f.22009.0.8"="genetic_principal_components_f22009_0_8",
    "f.22009.0.9"="genetic_principal_components_f22009_0_9",
    "f.22009.0.10"="genetic_principal_components_f22009_0_10",
    "f.22009.0.11"="genetic_principal_components_f22009_0_11",
    "f.22009.0.12"="genetic_principal_components_f22009_0_12",
    "f.22009.0.13"="genetic_principal_components_f22009_0_13",
    "f.22009.0.14"="genetic_principal_components_f22009_0_14",
    "f.22009.0.15"="genetic_principal_components_f22009_0_15",
    "f.22009.0.16"="genetic_principal_components_f22009_0_16",
    "f.22009.0.17"="genetic_principal_components_f22009_0_17",
    "f.22009.0.18"="genetic_principal_components_f22009_0_18",
    "f.22009.0.19"="genetic_principal_components_f22009_0_19",
    "f.22009.0.20"="genetic_principal_components_f22009_0_20"
)
names(phe) <- field_hash[names(phe)]

# The fields that we need for preprocessing (no PCs this time)
fields <- c(
    "eid",
    "sex_f31_0_0",
    "standing_height_f50_0_0",
    "standing_height_f50_1_0",
    "standing_height_f50_2_0",
    "standing_height_f50_3_0",
    "date_of_attending_assessment_centre_f53_0_0",
    "date_of_attending_assessment_centre_f53_1_0",
    "date_of_attending_assessment_centre_f53_2_0",
    "date_of_attending_assessment_centre_f53_3_0",
    "ethnic_background_f21000_0_0",
    "weight_f21002_0_0",
    "weight_f21002_1_0",
    "weight_f21002_2_0",
    "weight_f21002_3_0",
    "age_when_attended_assessment_centre_f21003_0_0",
    "genotype_measurement_batch_f22000_0_0",
    "genetic_sex_f22001_0_0"
)
phe <- phe[,fields,drop=FALSE]

# Define the number of samples for weight change case
dwc <- phe[,c(1,12:15)]
# Which visit has the most recorded cases?
nas <- apply(dwc,2,function(x) length(which(!is.na(x))))
#              eid weight_f21002_0_0 weight_f21002_1_0 weight_f21002_2_0
#           502510            499736             20303             42097
#weight_f21002_3_0
#              826
## Is seems that visit 2 is the most prominent
dwc13 <- dwc[,c(1,2,4)]
dwc13c <- dwc13[complete.cases(dwc13),]
## 42035 cases, not bad!

# Define the data to be further filtered
rownames(dwc13c) <- dwc13c$eid
rownames(phe) <- phe$eid
phe <- phe[rownames(dwc13c),,drop=FALSE]
# Still some sex missing
tmp <- phe[,c(2,11,16,17,18)]
cstmp <- complete.cases(tmp)
phe <- phe[cstmp,] # 40997 samples

# Primary filters in phenos:
# Keep:
# phe$sex_f31_0_0 == phe$genetic_sex_f22001_0_0
# phe$ethnic_background_f21000_0_0==1001
filt_gen <- phe$ethnic_background_f21000_0_0==1001 & 
    phe$sex_f31_0_0 == phe$genetic_sex_f22001_0_0
# dWC outliers?
diwc <- phe$weight_f21002_2_0 - phe$weight_f21002_0_0
m_diwc <- median(diwc)
i_diwc <- IQR(diwc)
diwc_filt <- diwc > m_diwc + 3*i_diwc | diwc < m_diwc - 3*i_diwc
# Filter
filt <- filt_gen & !diwc_filt # 37055 samples

# Final and also assign age^2, drop PCs, we won't need them
phe$delta_wc = diwc
#phe$age_when_attended_assessment_centre_f21003_0_0_sq <- 
#    phe$age_when_attended_assessment_centre_f21003_0_0^2
# For delta WC we also need height
#phe_keep <- phe[filt,-c(4:15,18)] 
#write.table(phe_keep,file="covariates_dwc.txt",quote=FALSE,sep="\t",
#    row.names=FALSE)
# For dWC betas we also need height
phe_keep_beta <- phe[filt,-c(4,6,8,10,11,13,15,18,19)]
#write.table(phe_keep_beta,file="covariates_beta_dwc.txt",quote=FALSE,sep="\t",
#    row.names=FALSE)

# Process kept data to match the format described by A. Eriksson:
data_processed <- phe_keep_beta %>%
  mutate(
    # Concatenate heights
    Height=paste(standing_height_f50_0_0,standing_height_f50_2_0,sep=";"),
    
    # Concatenate weights
    Weight=paste(weight_f21002_0_0,weight_f21002_2_0,sep=";"),
    
    # Calculate age at second assessment
    # Convert dates to Date format
    date_first=as.Date(date_of_attending_assessment_centre_f53_0_0),
    date_second=as.Date(date_of_attending_assessment_centre_f53_2_0),
    
    # Calculate difference in years
    years_diff=as.numeric(difftime(date_second,date_first,units="days"))/365.25,
    
    # Calculate age at second assessment
    age_second=age_when_attended_assessment_centre_f21003_0_0 + years_diff,
    
    # Concatenate initial age and new age
    Age=paste(age_when_attended_assessment_centre_f21003_0_0,round(age_second),
        sep=";"),
    Sex=sex_f31_0_0
  )

data_export <- data_processed[,c("eid","Sex","Age","Height","Weight")]
write.table(data_export,"covariates_beta_dwc.txt",sep="\t",row.names=FALSE,
    quote=FALSE)
```

### Preparation of the regression response for δWC

As the reponse variable when constructing the δWC PRS from summary statistics
was consisted of the Weight change slopes (betas) instead of actual delta 
values, we below perform the slope estimation for the UKB phenotypes prepared 
above following the process and scripts distributed by A. Eriksson:

```R
# Read in the prepared file
d  <- read.delim("covariates_beta_dwc.txt")
nr <- nrow(d)

# Init output variables
d$N <- d$M_Height <- d$M_Age <- d$S_Age <- d$M_Weight <- d$S_Weight <- 
    d$Weight_Alpha <- d$Weight_Beta <- rep(0,nr)

for (k in 1:nr) {
    # Split Age, BMI and Weight strings into vectors
    age_t <- suppressWarnings(as.numeric(strsplit(d$Age[k],";")[[1]]))
    weight <- suppressWarnings(as.numeric(strsplit(d$Weight[k],";")[[1]]))
    
    # Calculate mean height from available values
    if (is.numeric(d$Height[k])) # One value only
        d$M_Height[k] <- d$Height[k]
    else {
        h <- suppressWarnings(as.numeric(strsplit(d$Height[k],";")[[1]]))
        d$M_Height[k] <- mean(h,na.rm=TRUE)
    }

    # Select time points with BMI information and age in the right interval
    use <- !is.na(age_t) & !is.na(weight) & age_t > 18 #& age_t <= 65

    # Subest data to eligible time points and offset so that Alpha is estimate 
    # of Weight at age 18
    age_t <- age_t[use] - 18
    weight <- weight[use]

    d$N[k] <- length(age_t)
    if (d$N[k] >= 2) {
        d$M_Age[k] <- mean(age_t);
        d$M_Weight[k] <- mean(weight);
        
        # Unbiased estimator of standard error of age in the sample
        d$S_Age[k] <- sqrt(sum((age_t-mean(age_t))^2)/(d$N[k]-1));
        # Unbiased estimator of standard error of Weight in the sample
        d$S_Weight[k] = sqrt(sum((weight-mean(weight))^2)/(d$N[k]-1)); 
        
        # Models
        f = summary(lm(weight ~ age_t))
        
        # Estimator of intercept at age 18 - Weight
        d$Weight_Alpha[k] <- f$coefficients[1,1]
        # Estimator of slope of BMI
        d$Weight_Beta[k] <- f$coefficients[2,1]
    }
}

# Save data to file for individuals with >= 2 measurement points
cols_to_keep <- c("eid","Sex","N","M_Height","M_Age","S_Age","M_Weight",
    "S_Weight","Weight_Alpha","Weight_Beta")
d_keep <- d[d$N>=2,cols_to_keep]
write.table(d_keep,"weight_change_summary_adults.txt",
    sep="\t",quote=FALSE,row.names=FALSE)
```

### Preparation of UKB data with δWC samples - all available samples

#### Conversion to PLINK 1.90 format

First, get a list of samples to be kept:

```bash
awk 'NR>1 {print $1" "$1}' covariates_beta_dwc.txt > dwc_keep.txt
```

and then perform the conversion while keeping the samples with weight difference
available at the same time. `--memory` and `--threads` are adjusted for a 
machine with 64 available threads and 512GB RAM so as to process 22 chromosomes
concurrently. Adjust according to your environment.

```bash
for CHR in `seq 1 22`
do
  plink2 \
    --bgen ukb_chr${CHR}.bgen ref-first \
    --sample ukb_chr${CHR}.sample \
    --keep dwc_keep.txt \
    --snps-only \
    --max-alleles 2 \
    --rm-dup force-first \
    --make-bed \
    --memory 22000 \
    --threads 3 \
    --out ukb_dwc_chr${CHR} &
done
```

Now merge chromosomes for easier calculations. The resulting sizes after 
filtering samples for weight change are not prohibitive:

```bash
for CHR in `seq 1 22`
do
  echo ukb_dwc_chr${CHR} >> merge.list
done

plink \
  --merge-list merge.list \
  --make-bed \
  --out ukb_dwc
```

#### Basic SNP and sample filtering

Perform filtering with PLINK 1.9 to remove rare variants as well as lower sample 
and call rates. The following variant filters are recommended:

1. Variant call rate: >98%
2. Minor Allele Frequency: >0.05
3. Hardy-Weinberg equilibrium p-value: >10<sup>-6</sup>

The following sample filters are recommended:

1. Sample call rate: > 95%
2. Heterozygosity: median(heterozygosity) &plusmn; 3 &times; IQR
3. Identity By Descent: > 0.5 (optional)

Below, we sequentially apply variant and sample filters according to widely
accepted [best practices](https://onlinelibrary.wiley.com/doi/10.1002/sim.6605).

First, we calculate heterozygosities:

```bash
plink \
  --bfile ukb_dwc \
  --out ukb_dwc --het
```

Then, we create a file with sample names with heterozygosities within the limits
of our filter (sample filter #2 above):

```r
Rscript \
  -e '{
    hetData <- read.table("ukb_dwc.het",header=TRUE,
        check.names=FALSE)
    rownames(hetData) <- hetData$IID
    heterozygosity <- 1 - hetData[,3]/hetData[,5]
    names(heterozygosity) <- rownames(hetData)
    avg <- median(heterozygosity,na.rm=TRUE)
    dev <- IQR(heterozygosity,na.rm=TRUE)
    keep <- heterozygosity > avg - 3*dev & heterozygosity < avg + 3*dev
    goodHet <- rownames(hetData)[keep]
    write.table(hetData[keep,c("FID","IID"),drop=FALSE],
      file="het_samples_pass.txt",col.names=FALSE,row.names=FALSE,quote=FALSE)
  }'
```

The file `het_samples_pass.txt` will be used after the next filters to compile 
a final list of samples to be kept in later analysis.

The following `plink` command will apply variant filters #1,2,3 and sample 
filter #1:

```bash
plink \
  --bfile ukb_dwc \
  --out ukb_dwc_tmp \
  --make-bed \
  --geno 0.02 \
  --maf 0.05 \
  --hwe 0.000001 \
  --mind 0.05
```

We now create files with variants and samples to *keep* (samples to keep are 
merged with those passing heterozygosity filters):

```r
cut -d" " -f1-2 ukb_dwc_tmp.fam > generic_samples_pass.txt
cut -f2 ukb_dwc_tmp.bim > generic_variants_pass.txt

Rscript \
  -e '{
    het <- read.table("het_samples_pass.txt")
    rownames(het) <- het[,2]
    gen <- read.table("generic_samples_pass.txt")
    rownames(gen) <- gen[,2]
    pass <- intersect(rownames(het),rownames(gen))
    write.table(gen[pass,,drop=FALSE],file="all_samples_pass.txt",
      col.names=FALSE,row.names=FALSE,quote=FALSE)
  }'
```

Based on the variant and sample content of the files `generic_variants_pass.txt`
and `all_samples_pass.txt` we create filtered PLINK files. These will be used
for kinship analysis filtering (can be skipped if not necesseary) and PCA.

```bash
plink \
  --bfile ukb_dwc \
  --out ukb_dwc_filtered \
  --extract generic_variants_pass.txt \
  --keep all_samples_pass.txt \
  --make-bed
```

Kinship analysis is performed with KING (no LD pruning):

```bash
king \
  -b ukb_dwc_filtered.bed \
  --unrelated \
  --degree 2 \
  --prefix individuals_
```

The aforementioned command will create the file `individuals_unrelated.txt` 
which we use with PLINK to create the dataset with related individuals removed
(optional).

```bash
plink \
  --bfile ukb_dwc_filtered \
  --keep individuals_unrelated.txt \
  --out ukb_dwc_final \
  --make-bed
```

#### Perform LD pruning only δBMI data

Prior to PCA, it is common practice to not use variants in LD, so LD pruning is
performed.

* `--indep-pairwise [window size] [step size/variant count)] [R2]` controls the 
pruning parameters. e.g. indep `50 5 0.9` and generates a list of markers in 
approximate linkage equilibrium - takes 50 SNPs at a time and then shifts by 5 
for the window. R2 is the cut-off for linkage disequilibrium.
* We follow an iterative procedure similar to the
[FinnGen](https://finngen.gitbook.io/documentation/methods/phewas/quality-checks)
consortium, where R2 is decreased by a constant step until ~200,000 variants are
included in the pruned set.

```bash
mkdir pruned

TOTAL=10000000
R2=0.95
STEP=0.05

while [ $TOTAL -gt 200000 ]
do
  R2=$(echo "$R2-$STEP" | bc)
  
  plink \
    --bfile ukb_dwc_final \
    --indep-pairwise 500 50 $R2 \
    --out ./pruned/ukb_dwc_final
    
  TOTAL=$(wc -l ./pruned/*.in | awk '{print $1}')
  
  echo "Parameters --indep-pairwise 500 50 $R2 yield $TOTAL SNPs" >> \
    pruning.log
done
```

The final value of `$R2` is `0.3` and the number of SNPs is 187,959, so we
perform the final pruning:

```bash
rm -r ./pruned/*

FINAL_R2=`tail -1 pruning.log | awk '{print $(NF-3)}' | sed 's/^\./0./'`
plink \
  --bfile ukb_dwc_final \
  --indep-pairwise 500 50 $FINAL_R2 \
  --out ukb_dwc_final_pruned

plink \
  --bfile ukb_dwc_final \
  --extract ukb_dwc_final_pruned.prune.in \
  --make-bed \
  --out ukb_dwc_final_pruned

rm *.nosex *.log
```

#### Principal Component Analysis with FlashPca2

Perform PCA using the FlashPCA2:

```
flashpca \
  --bfile ukb_dwc_final_pruned \
  --ndim 20 \
  --outpc pcs_ukb_dwc.txt \
  --outvec eig_ukb_dwc.txt \
  --outload loads_ukb_dwc.txt
```

#### Merge PCs with covariates file

The process is completed in R:

```R
covars <- read.delim("weight_change_summary_adults.txt")
rownames(covars) <- covars$eid
covars <- covars[,-c(3,6:9)]
names(covars)[3:4] <- c("Height","Age")
pcs <- read.delim("pcs_ukb_dwc.txt")
rownames(pcs) <- pcs$IID
covars <- covars[rownames(pcs),,drop=FALSE]
covars20 <- cbind(covars,pcs[,3:ncol(pcs)])
covars10 <- cbind(covars,pcs[,3:12])
write.table(covars20,file="covariates_dwc_20pc.txt",sep="\t",quote=FALSE,
    row.names=FALSE)
write.table(covars10,file="covariates_dwc_10pc.txt",sep="\t",quote=FALSE,
    row.names=FALSE)
```

### Assessment of Polygenic Risk Scores for δΒΜΙ

We assume the same workspace as before.

#### All-inclusive cohorts

As previously [described](https://github.com/moulos-lab/better4u/tree/main/03_prs_derivation#weight-change), we have derived 5 versions of a PRS using all available participating cohorts. The 
files are:

- `b4u_wc_prscs_original.prs`
- `b4u_wc_prscs_robust.prs`
- `b4u_wc_prscs_recalibrated.prs`
- `b4u_wc_sbrc_tgp.prs`
- `b4u_wc_sbrc_ukb.prs`

We can retrieve these files from Google Drive:

```bash
# b4u_wc_prscs_original.prs
gdown 15136PvuqV9hpKV-Bc5LZTGv_FUv_Ew1H --output b4u_wc_prscs_original.prs

# b4u_wc_prscs_robust.prs
gdown 1c2OmdJJKtgvG8-gzwSHZLUHc0A0ad-hW --output b4u_wc_prscs_robust.prs

# b4u_wc_prscs_recalibrated.prs
gdown 17rrvwqluPgnl_fQPeIgFNsuwkFs9SaHB --output b4u_wc_prscs_recalibrated.prs

# b4u_wc_sbrc_tgp.prs
gdown 1LY1OZOdSCM3eyfMGxeQpK2JW1j40gAv1 --output b4u_wc_sbrc_tgp.prs

# b4u_wc_sbrc_ukb.prs
gdown 1hFGxhncwmY4lTMAOtS7wa5Hx7cOyngsC --output b4u_wc_sbrc_ukb.prs
``` 

Below we produce evaluation metrics for UKB after basic phenotypic QC. To 
achieve maximum PRS coverage we skip QC in genotypic data.

```R
WORKSPACE <- Sys.getenv("WORKSPACE")

source(file.path(WORKSPACE,"better4u","03_prs_derivation","evalfuns.R"))

trait <- "Weight_Beta"
genoBase <- "ukb_dwc"
covFile <- "covariates_dwc_20pcs.txt"
#covFile <- "covariates_dwc_10pcs.txt"

# PRS-CS with 1000 genomes LD - original version
prsFile_PRSCS_ORG <- "b4u_wc_prscs_original.prs"
# We already have BIMs from the preprocessing
sanFile_PRSCS_ORG <- sanitizePrs(prsFile_PRSCS_ORG,genoBase,perChr=TRUE,
    from="ready",rc=0.4)
M_PRSCS_ORG <- evalPrs(sanFile_PRSCS_ORG,covFile,trait,genoBase,iidCol=1,
    perChr=TRUE,remFile="removed.txt",rc=0.4)

# PRS-CS with 1000 genomes LD - robust SNPs - original weights
prsFile_PRSCS_ROB <- "b4u_wc_prscs_robust.prs"
sanFile_PRSCS_ROB <- sanitizePrs(prsFile_PRSCS_ROB,genoBase,perChr=TRUE,
    from="ready",rc=0.4)
M_PRSCS_ROB <- evalPrs(sanFile_PRSCS_ROB,covFile,trait,genoBase,iidCol=1,
    perChr=TRUE,remFile="removed.txt",rc=0.4)

# PRS-CS with 1000 genomes LD - robust SNPs - original weights
prsFile_PRSCS_REC <- "b4u_wc_prscs_recalibrated.prs"
sanFile_PRSCS_REC <- sanitizePrs(prsFile_PRSCS_REC,genoBase,perChr=TRUE,
    from="ready",rc=0.4)
M_PRSCS_REC <- evalPrs(sanFile_PRSCS_REC,covFile,trait,genoBase,iidCol=1,
    perChr=TRUE,remFile="removed.txt",rc=0.4)

# SBayesRC with custom 1000 genomes LD
prsFile_SBRC_TGP <- "b4u_wc_sbrc_tgp.prs"
sanFile_SBRC_TGP <- sanitizePrs(prsFile_SBRC_TGP,genoBase,perChr=TRUE,
    from="ready",rc=0.4)
M_SBRC_TGP <- evalPrs(sanFile_SBRC_TGP,covFile,trait,genoBase,iidCol=1,
    perChr=TRUE,remFile="removed.txt",rc=0.4)

# SBayesRC with built-in UKB LD
prsFile_SBRC_UKB <- "b4u_wc_sbrc_ukb.prs"
sanFile_SBRC_UKB <- sanitizePrs(prsFile_SBRC_UKB,genoBase,perChr=TRUE,
    from="ready",rc=0.4)
M_SBRC_UKB <- evalPrs(sanFile_SBRC_UKB,covFile,trait,genoBase,iidCol=1,
    perChr=TRUE,remFile="removed.txt",rc=0.4)
    
save(M_PRSCS_ORG,M_PRSCS_ROB,M_PRSCS_REC,M_SBRC_TGP,M_SBRC_UKB,
    file="metrics_nopc_onlybeta.rda")
#save(M_PRSCS_ORG,M_PRSCS_ROB,M_PRSCS_REC,M_SBRC_TGP,M_SBRC_UKB,
#    file="metrics_10pc.rda")

# Some additional values regarding coverage
N_PRSCS_ORG <- as.numeric(countLines(prsFile_PRSCS_ORG)) - 1
N_PRSCS_ROB <- as.numeric(countLines(prsFile_PRSCS_ROB)) - 1
N_PRSCS_REC <- as.numeric(countLines(prsFile_PRSCS_REC)) - 1
N_SBRC_TGP <- as.numeric(countLines(prsFile_SBRC_TGP)) - 1
N_SBRC_UKB <- as.numeric(countLines(prsFile_SBRC_UKB)) - 1

# Gather metrics to share and PRS (just values, no ids of any kind)
metrics <- data.frame(
    original=formatC(M_PRSCS_ORG$metrics,digits=6),
    robust=formatC(M_PRSCS_ROB$metrics,digits=6),
    recalibrated=formatC(M_PRSCS_REC$metrics,digits=6),
    sbayes_tgp=formatC(M_SBRC_TGP$metrics,digits=6),
    sbayes_ukb=formatC(M_SBRC_UKB$metrics,digits=6)
)

# Additonal data
add <- as.data.frame(rbind(
    as.integer(c(N_PRSCS_ORG,N_PRSCS_ROB,N_PRSCS_REC,N_SBRC_TGP,N_SBRC_UKB)),
    round(100*as.numeric(metrics[nrow(metrics),])/
        c(N_PRSCS_ORG,N_PRSCS_ROB,N_PRSCS_REC,N_SBRC_TGP,N_SBRC_UKB),digits=2)
))
rownames(add) <- c("snps_total","coverage")
colnames(add) <- names(metrics)

# Final metrics data frame
finalMetrics <- rbind(metrics,add)

write.table(finalMetrics,file="b4u_wc_prs_metrics_UKB_20pc.txt",sep="\t",
    quote=FALSE,col.names=NA)
#write.table(finalMetrics,file="b4u_wc_prs_metrics_UKB_10pc.txt",sep="\t",
#    quote=FALSE,col.names=NA)
```

#### LOO cohorts

*WIP*

### Preparation of UKB data with δWC samples - 18-65 years

#### Conversion to PLINK 1.90 format

First, get a list of sample to be kept:

```bash
awk 'NR>1 {print $1" "$1}' weight_summary_adults_18_65.txt > dwc1865_keep.txt
```

and then perform the conversion while keeping the samples with weight difference
available at the same time. `--memory` and `--threads` are adjusted for a 
machine with 64 available threads and 512GB RAM so as to process 22 chromosomes
concurrently. Adjust according to your environment.

```bash
for CHR in `seq 1 22`
do
  plink2 \
    --bgen ukb_chr${CHR}.bgen ref-first \
    --sample ukb_chr${CHR}.sample \
    --keep dwc1865_keep.txt \
    --snps-only \
    --max-alleles 2 \
    --rm-dup force-first \
    --make-bed \
    --memory 22000 \
    --threads 3 \
    --out ukb_dwc1865_chr${CHR} &
done
```

Now merge chromosomes for easier calculations. The resulting sizes after 
filtering samples for weight change are not prohibitive:

```bash
for CHR in `seq 1 22`
do
  echo ukb_dwc1865_chr${CHR} >> merge.list
done

plink \
  --merge-list merge.list \
  --make-bed \
  --out ukb_dwc1865
```

#### Basic SNP and sample filtering

Perform filtering with PLINK 1.9 to remove rare variants as well as lower sample 
and call rates. The following variant filters are recommended:

1. Variant call rate: >98%
2. Minor Allele Frequency: >0.05
3. Hardy-Weinberg equilibrium p-value: >10<sup>-6</sup>

The following sample filters are recommended:

1. Sample call rate: > 95%
2. Heterozygosity: median(heterozygosity) &plusmn; 3 &times; IQR
3. Identity By Descent: > 0.5 (optional)

Below, we sequentially apply variant and sample filters according to widely
accepted [best practices](https://onlinelibrary.wiley.com/doi/10.1002/sim.6605).

First, we calculate heterozygosities:

```bash
plink \
  --bfile ukb_dwc1865 \
  --out ukb_dwc1865 --het
```

Then, we create a file with sample names with heterozygosities within the limits
of our filter (sample filter #2 above):

```r
Rscript \
  -e '{
    hetData <- read.table("ukb_dwc1865.het",header=TRUE,
        check.names=FALSE)
    rownames(hetData) <- hetData$IID
    heterozygosity <- 1 - hetData[,3]/hetData[,5]
    names(heterozygosity) <- rownames(hetData)
    avg <- median(heterozygosity,na.rm=TRUE)
    dev <- IQR(heterozygosity,na.rm=TRUE)
    keep <- heterozygosity > avg - 3*dev & heterozygosity < avg + 3*dev
    goodHet <- rownames(hetData)[keep]
    write.table(hetData[keep,c("FID","IID"),drop=FALSE],
      file="het_samples_pass.txt",col.names=FALSE,row.names=FALSE,quote=FALSE)
  }'
```

The file `het_samples_pass.txt` will be used after the next filters to compile 
a final list of samples to be kept in later analysis.

The following `plink` command will apply variant filters #1,2,3 and sample 
filter #1:

```bash
plink \
  --bfile ukb_dwc1865 \
  --out ukb_dwc1865_tmp \
  --make-bed \
  --geno 0.02 \
  --maf 0.05 \
  --hwe 0.000001 \
  --mind 0.05
```

We now create files with variants and samples to *keep* (samples to keep are 
merged with those passing heterozygosity filters):

```r
cut -d" " -f1-2 ukb_dwc1865_tmp.fam > generic_samples_pass.txt
cut -f2 ukb_dwc1865_tmp.bim > generic_variants_pass.txt

Rscript \
  -e '{
    het <- read.table("het_samples_pass.txt")
    rownames(het) <- het[,2]
    gen <- read.table("generic_samples_pass.txt")
    rownames(gen) <- gen[,2]
    pass <- intersect(rownames(het),rownames(gen))
    write.table(gen[pass,,drop=FALSE],file="all_samples_pass.txt",
      col.names=FALSE,row.names=FALSE,quote=FALSE)
  }'
```

Based on the variant and sample content of the files `generic_variants_pass.txt`
and `all_samples_pass.txt` we create filtered PLINK files. These will be used
for kinship analysis filtering (can be skipped if not necesseary) and PCA.

```bash
plink \
  --bfile ukb_dwc1865 \
  --out ukb_dwc1865_filtered \
  --extract generic_variants_pass.txt \
  --keep all_samples_pass.txt \
  --make-bed
```

Kinship analysis is performed with KING (no LD pruning):

```bash
king \
  -b ukb_dwc1865_filtered.bed \
  --unrelated \
  --degree 2 \
  --prefix individuals_
```

The aforementioned command will create the file `individuals_unrelated.txt` 
which we use with PLINK to create the dataset with related individuals removed
(optional).

```bash
plink \
  --bfile ukb_dwc1865_filtered \
  --keep individuals_unrelated.txt \
  --out ukb_dwc1865_final \
  --make-bed
```

#### Perform LD pruning only δBMI data

Prior to PCA, it is common practice to not use variants in LD, so LD pruning is
performed.

* `--indep-pairwise [window size] [step size/variant count)] [R2]` controls the 
pruning parameters. e.g. indep `50 5 0.9` and generates a list of markers in 
approximate linkage equilibrium - takes 50 SNPs at a time and then shifts by 5 
for the window. R2 is the cut-off for linkage disequilibrium.
* We follow an iterative procedure similar to the
[FinnGen](https://finngen.gitbook.io/documentation/methods/phewas/quality-checks)
consortium, where R2 is decreased by a constant step until ~200,000 variants are
included in the pruned set.

```bash
mkdir pruned

TOTAL=10000000
R2=0.95
STEP=0.05

while [ $TOTAL -gt 200000 ]
do
  R2=$(echo "$R2-$STEP" | bc)
  
  plink \
    --bfile ukb_dwc_final \
    --indep-pairwise 500 50 $R2 \
    --out ./pruned/ukb_dwc_final
    
  TOTAL=$(wc -l ./pruned/*.in | awk '{print $1}')
  
  echo "Parameters --indep-pairwise 500 50 $R2 yield $TOTAL SNPs" >> \
    pruning.log
done
```

The final value of `$R2` is `0.3` and the number of SNPs is 187,959, so we
perform the final pruning:

```bash
rm -r ./pruned/*

FINAL_R2=`tail -1 pruning.log | awk '{print $(NF-3)}' | sed 's/^\./0./'`
plink \
  --bfile ukb_dwc1865_final \
  --indep-pairwise 500 50 $FINAL_R2 \
  --out ukb_dwc1865_final_pruned

plink \
  --bfile ukb_dwc1865_final \
  --extract ukb_dwc1865_final_pruned.prune.in \
  --make-bed \
  --out ukb_dwc1865_final_pruned

rm *.nosex *.log
```

#### Principal Component Analysis with FlashPca2

Perform PCA using the FlashPCA2:

```
flashpca \
  --bfile ukb_dwc1865_final_pruned \
  --ndim 20 \
  --outpc pcs_ukb_dwc1865.txt \
  --outvec eig_ukb_dwc1865.txt \
  --outload loads_ukb_dwc1865.txt
```

#### Merge PCs with covariates file

The process is completed in R:

```R
covars <- read.delim("weight_summary_adults_18_65.txt")
rownames(covars) <- covars$eid
covars <- covars[,-c(3:6,10:15)]
pcs <- read.delim("pcs_ukb_dwc1865.txt")
rownames(pcs) <- pcs$IID
covars <- covars[rownames(pcs),,drop=FALSE]
covars20 <- cbind(covars,pcs[,3:ncol(pcs)])
covars10 <- cbind(covars,pcs[,3:12])
write.table(covars20,file="covariates_dwc1865_20pcs.txt",sep="\t",quote=FALSE,
    row.names=FALSE)
write.table(covars10,file="covariates_dwc1865_10pcs.txt",sep="\t",quote=FALSE,
    row.names=FALSE)
```

### Assessment of Polygenic Risk Scores for weight change

We assume the same workspace as before.

#### All-inclusive cohorts

As previously [described](https://github.com/moulos-lab/better4u/tree/main/03_prs_derivation#weight-change), we have derived 5 versions of a PRS using all available participating cohorts. The 
files are:

- `b4u_wc_prscs_original.prs`
- `b4u_wc_prscs_robust.prs`
- `b4u_wc_prscs_recalibrated.prs`
- `b4u_wc_sbrc_tgp.prs`
- `b4u_wc_sbrc_ukb.prs`

We can retrieve these files from Google Drive:

```bash
# b4u_wc_prscs_original.prs
gdown 15136PvuqV9hpKV-Bc5LZTGv_FUv_Ew1H --output b4u_wc_prscs_original.prs

# b4u_wc_prscs_robust.prs
gdown 1c2OmdJJKtgvG8-gzwSHZLUHc0A0ad-hW --output b4u_wc_prscs_robust.prs

# b4u_wc_prscs_recalibrated.prs
gdown 17rrvwqluPgnl_fQPeIgFNsuwkFs9SaHB --output b4u_wc_prscs_recalibrated.prs

# b4u_wc_sbrc_tgp.prs
gdown 1LY1OZOdSCM3eyfMGxeQpK2JW1j40gAv1 --output b4u_wc_sbrc_tgp.prs

# b4u_wc_sbrc_ukb.prs
gdown 1hFGxhncwmY4lTMAOtS7wa5Hx7cOyngsC --output b4u_wc_sbrc_ukb.prs
``` 

Below we produce evaluation metrics for UKB after basic phenotypic QC. To 
achieve maximum PRS coverage we skip QC in genotypic data.

```R
WORKSPACE <- Sys.getenv("WORKSPACE")

source(file.path(WORKSPACE,"better4u","03_prs_derivation","evalfuns.R"))

trait <- "Weight_Beta"
genoBase <- "ukb_dwc1865"
covFile <- "covariates_dwc1865_20pcs.txt"
#covFile <- "covariates_dwc1865_10pcs.txt"

# PRS-CS with 1000 genomes LD - original version
prsFile_PRSCS_ORG <- "b4u_wc_prscs_original.prs"
# We already have BIMs from the preprocessing
sanFile_PRSCS_ORG <- sanitizePrs(prsFile_PRSCS_ORG,genoBase,perChr=TRUE,
    from="ready",rc=0.4)
M_PRSCS_ORG <- evalPrs(sanFile_PRSCS_ORG,covFile,trait,genoBase,iidCol=1,
    perChr=TRUE,remFile="removed.txt",rc=0.4)

# PRS-CS with 1000 genomes LD - robust SNPs - original weights
prsFile_PRSCS_ROB <- "b4u_wc_prscs_robust.prs"
# We have the BIMs, no need for bgen or ukb to be TRUE
sanFile_PRSCS_ROB <- sanitizePrs(prsFile_PRSCS_ROB,genoBase,perChr=TRUE,
    from="ready",rc=0.4)
M_PRSCS_ROB <- evalPrs(sanFile_PRSCS_ROB,covFile,trait,genoBase,iidCol=1,
    perChr=TRUE,remFile="removed.txt",rc=0.4)

# PRS-CS with 1000 genomes LD - robust SNPs - original weights
prsFile_PRSCS_REC <- "b4u_wc_prscs_recalibrated.prs"
sanFile_PRSCS_REC <- sanitizePrs(prsFile_PRSCS_REC,genoBase,perChr=TRUE,
    from="ready",rc=0.4)
M_PRSCS_REC <- evalPrs(sanFile_PRSCS_REC,covFile,trait,genoBase,iidCol=1,
    perChr=TRUE,remFile="removed.txt",rc=0.4)

# SBayesRC with custom 1000 genomes LD
prsFile_SBRC_TGP <- "b4u_wc_sbrc_tgp.prs"
sanFile_SBRC_TGP <- sanitizePrs(prsFile_SBRC_TGP,genoBase,perChr=TRUE,
    from="ready",rc=0.4)
M_SBRC_TGP <- evalPrs(sanFile_SBRC_TGP,covFile,trait,genoBase,iidCol=1,
    perChr=TRUE,remFile="removed.txt",rc=0.4)

# SBayesRC with built-in UKB LD
prsFile_SBRC_UKB <- "b4u_wc_sbrc_ukb.prs"
sanFile_SBRC_UKB <- sanitizePrs(prsFile_SBRC_UKB,genoBase,perChr=TRUE,
    from="ready",rc=0.4)
M_SBRC_UKB <- evalPrs(sanFile_SBRC_UKB,covFile,trait,genoBase,iidCol=1,
    perChr=TRUE,remFile="removed.txt",rc=0.4)
    
save(M_PRSCS_ORG,M_PRSCS_ROB,M_PRSCS_REC,M_SBRC_TGP,M_SBRC_UKB,
    file="metrics_nopc_onlybeta.rda")
#save(M_PRSCS_ORG,M_PRSCS_ROB,M_PRSCS_REC,M_SBRC_TGP,M_SBRC_UKB,
#    file="metrics_10pc.rda")

# Some additional values regarding coverage
N_PRSCS_ORG <- as.numeric(countLines(prsFile_PRSCS_ORG)) - 1
N_PRSCS_ROB <- as.numeric(countLines(prsFile_PRSCS_ROB)) - 1
N_PRSCS_REC <- as.numeric(countLines(prsFile_PRSCS_REC)) - 1
N_SBRC_TGP <- as.numeric(countLines(prsFile_SBRC_TGP)) - 1
N_SBRC_UKB <- as.numeric(countLines(prsFile_SBRC_UKB)) - 1

# Gather metrics to share and PRS (just values, no ids of any kind)
metrics <- data.frame(
    original=formatC(M_PRSCS_ORG$metrics,digits=6),
    robust=formatC(M_PRSCS_ROB$metrics,digits=6),
    recalibrated=formatC(M_PRSCS_REC$metrics,digits=6),
    sbayes_tgp=formatC(M_SBRC_TGP$metrics,digits=6),
    sbayes_ukb=formatC(M_SBRC_UKB$metrics,digits=6)
)

# Additonal data
add <- as.data.frame(rbind(
    as.integer(c(N_PRSCS_ORG,N_PRSCS_ROB,N_PRSCS_REC,N_SBRC_TGP,N_SBRC_UKB)),
    round(100*as.numeric(metrics[nrow(metrics),])/
        c(N_PRSCS_ORG,N_PRSCS_ROB,N_PRSCS_REC,N_SBRC_TGP,N_SBRC_UKB),digits=2)
))
rownames(add) <- c("snps_total","coverage")
colnames(add) <- names(metrics)

# Final metrics data frame
finalMetrics <- rbind(metrics,add)

write.table(finalMetrics,file="b4u_wc_prs_metrics_UKB_nopc_onlybeta.txt",sep="\t",
    quote=FALSE,col.names=NA)
#write.table(finalMetrics,file="b4u_wc_prs_metrics_UKB_10pc.txt",sep="\t",
#    quote=FALSE,col.names=NA)
```

#### LOO cohorts

*WIP*
