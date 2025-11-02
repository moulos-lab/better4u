BETTER4U Task 3.7 script distribution and data collection
================================================================================

## Authors

**Panagiotis Moulos** <sup>1,2</sup><br/>
<sup>1</sup>Harokopio University of Athens, <sup>2</sup>BSRC Alexander Fleming

### Contact details

* Panagiotis Moulos (moulos@fleming.gr)
* Jon Anders Eriksson (anders.eriksson@ut.ee)

# Introduction

This document explains the data that must be collected by partners for the 
completion of BETTER4U WP3 Task 3.7: Development of an iterative PRS extraction 
pipeline with the purpose of collecting data for a potential publication, apart
from performance metrics.

The data that will be collected from each cohort can be summarized to:

* Calculation of baseline PRS for BMI using PRS-CS and SBayesRC derived scores
* Calculation of baseline PRS for weight change using PRS-CS(x) derived score (_WIP_)
* Calculation of baseline PRS for BMI using robustified set of SNPs from PRS-CS bootstrap
* Calculation of performance metrics for all PRS versions

For the completion of the above tasks, a set of files will be distributed along
with instructions of how to use them. Please make sure that you maintain your
cohort data also as PLINK BED files, either as a whole or per chromosome.

## Required software

* [PLINK 1.90](https://www.cog-genomics.org/plink/)
* [R](https://www.r-project.org/)
* [gdown](https://github.com/wkentaro/gdown) (optional)
* [This](https://github.com/moulos-lab/better4u/blob/main/03_prs_derivation/evalfuns.R) set of R scripts

## Requested work description

In the following, we assume that the prefix name of your cohort is `COHORT`. For
the `COHORT` we suppose that one of the following applies:

a. There is one triplet (BED, BIM, FAM) of PLINK files for all the variants in
your cohort named:

```
COHORT.bed
COHORT.bim
COHORT.fam
```

b. There is one triplet of PLINK files per chromosome, named:

```
COHORT_chrZ.bed
COHORT_chrZ.bim
COHORT_chrZ.fam
```

where `Z = 1..22`. **Note** that the separator between `COHORT` and `chrZ`
should be an underscore (`_`). Altough in the accompanying scripts the separator
can be specified, we suggest to keep it to undescore (`_`). You can explore the
`sanitizePrs()` or `evalPrs()` functions in `evalfuns.R`.

### Calculation of baseline PRS for BMI and weight change

1. Setup the working folder. This can be a folder of your preference in your
system, e.g. `/home/myself/better4u/prs`:

```bash
# The prepending export is important, please do not omit!
export WORKSPACE=/home/myself/better4u/prs
mkdir -p $WORKSPACE && cd $WORKSPACE

git clone https://github.com/moulos-lab/better4u.git
```

2. Download the 5 PRS score files from [here](https://drive.google.com/drive/folders/14mwfzyHly-FSReZbpFojhsiNT6x2ph-f?usp=sharing). The folder contains the following files:

- `b4u_bmi_prscs_original.prs`: PRS constructed with PRS-CS with built-in EUR 
1000 genomes panel and a BIM file constructed from the BETTER4U BMI summary
statistics
- `b4u_bmi_prscs_robust.prs`: PRS-CS with built-in EUR 1000 genomes panel and 
robust SNP set after bootstrap and original PRS-CS betas
- `b4u_bmi_prscs_recalibrated.prs`: PRS-CS with built-in EUR 1000 genomes panel,
robust SNP set and recalibrated betas with PRS-CS
- `b4u_bmi_sbrc_tgp.prs`: GCTB SBayesRC with custom LD panel from 1000 genomes 
EUR samples
- `b4u_bmi_sbrc_ukb.prs`: GCTB SBayesRC with built-in UKB LD panel

If you have [gdown](https://github.com/wkentaro/gdown) installed:

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

Otherwise, just download from Google Drive and copy to the the desired location.

3. Calculate baseline PRSs based on the file from (1). We assume that all the
covariate files used in fastGWA analysis are placed in the file `covariates.txt` 
and is present in `$WORKSPACE`. The same covariates as the ones used with GCTA 
`--fastGWA` should be used. Note that the trait should also be attached to this
`covariates.txt` file. Then, within R:

```R
WORKSPACE <- Sys.getenv("WORKSPACE")

source(file.path(WORKSPACE,"better4u","03_prs_derivation","evalfuns.R"))

# Genetic data file (whether single or per chromosome)
genoBase <- "COHORT"

# Covariates file (sex, age, age^2, PCs) and trait
covFile <- "covariates.txt"

# Assuming this is the BMI name in covariates, otherwise adjust accordingly
trait <- "bmi" 

# Then, with the downloaded PRS score files

# PRS-CS with built-in EUR 1000 genomes panel
prsFile_PRSCS_ORG <- "b4u_bmi_prscs_original.prs"

# PRS-CS with built-in EUR 1000 genomes panel - robust SNPs original betas
prsFile_PRSCS_ROB <- "b4u_bmi_prscs_robust.prs"

# PRS-CS with built-in EUR 1000 genomes panel - recalibrated
prsFile_PRSCS_REC <- "b4u_bmi_prscs_recalibrated.prs"

# GCTB SBayesRC with TGP panel
prsFile_SBRC_TGP <- "b4u_bmi_sbrc_tgp.prs"

# GCTB SBayesRC with built-in LD panel
prsFile_SBRC_UKB <- "b4u_bmi_sbrc_ukb.prs"

# Sanitize the PRS scores (allele flip checking, coverage, etc)
# Single PLINK file:
sanFile_PRSCS_ORG <- sanitizePrs(prsFile_PRSCS_ORG,genoBase,from="ready")
sanFile_PRSCS_ROB <- sanitizePrs(prsFile_PRSCS_ROB,genoBase,from="ready")
sanFile_PRSCS_REC <- sanitizePrs(prsFile_PRSCS_REC,genoBase,from="ready")
sanFile_SBRC_TGP <- sanitizePrs(prsFile_SBRC_TGP,genoBase,from="sbayesrc")
sanFile_SBRC_UKB <- sanitizePrs(prsFile_SBRC_UKB,genoBase,from="sbayesrc")

# PLINK files per chromosome - if required
#sanFile_PRSCS_ORG <- sanitizePrs(prsFile_PRSCS_ORG,genoBase,perChr=TRUE,
#   from="ready",rc=0.2)
#sanFile_PRSCS_ROB <- sanitizePrs(prsFile_PRSCS_ROB,genoBase,perChr=TRUE,
#   from="ready",rc=0.2)
#sanFile_PRSCS_REC <- sanitizePrs(prsFile_PRSCS_REC,genoBase,perChr=TRUE,
#   from="ready",rc=0.2)
#sanFile_SBRC_TGP <- sanitizePrs(prsFile_SBRC_TGP,genoBase,perChr=TRUE,
#   from="sbayesrc",rc=0.2)
#sanFile_SBRC_UKB <- sanitizePrs(prsFile_SBRC_UKB,genoBase,perChr=TRUE,
#   from="sbayesrc",rc=0.2)

# Five files are written:
# b4u_bmi_prscs_original.prs.san
# b4u_bmi_prscs_robust.prs.san
# b4u_bmi_prscs_recalibrated.prs.san
# b4u_bmi_sbrc_tgp.prs.san
# b4u_bmi_sbrc_ukb.prs.san

# Then feed to the PRS/metrics function calculation
# Single PLINK file:
M_PRSCS_ORG <- evalPrs(sanFile_PRSCS_ORG,covFile,trait,genoBase)
M_PRSCS_ROB <- evalPrs(sanFile_PRSCS_ROB,covFile,trait,genoBase)
M_PRSCS_REC <- evalPrs(sanFile_PRSCS_REC,covFile,trait,genoBase)
M_SBRC_TGP <- evalPrs(sanFile_SBRC_TGP,covFile,trait,genoBase)
M_SBRC_UKB <- evalPrs(sanFile_SBRC_UKB,covFile,trait,genoBase)

# PLINK files per chromosome - if required
#M_PRSCS_ORG <- evalPrs(sanFile_PRSCS_ORG,covFile,trait,genoBase,
#   perChr=TRUE,rc=0.2)
#M_PRSCS_ROB <- evalPrs(sanFile_PRSCS_ROB,covFile,trait,genoBase,
#   perChr=TRUE,rc=0.2)
#M_PRSCS_REC <- evalPrs(sanFile_PRSCS_REC,covFile,trait,genoBase,
#   perChr=TRUE,rc=0.2)
#M_SBRC_TGP <- evalPrs(sanFile_SBRC_TGP,covFile,trait,genoBase,
#   perChr=TRUE,rc=0.2)
#M_SBRC_UKB <- evalPrs(sanFile_SBRC_UKB,covFile,trait,genoBase,
#   perChr=TRUE,rc=0.2)

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

write.table(finalMetrics,file="b4u_bmi_prs_metrics.txt",sep="\t",
    quote=FALSE,col.names=NA)

# Now write the actual PRS values for each version
prs_PRSCS_ORG <- as.data.frame(round(M_PRSCS_ORG$prs,3))
prs_PRSCS_ROB <- as.data.frame(round(M_PRSCS_ROB$prs,3))
prs_PRSCS_REC <- as.data.frame(round(M_PRSCS_REC$prs,3))
prs_SBRC_TGP <- as.data.frame(round(M_SBRC_TGP$prs,3))
prs_SBRC_UKB <- as.data.frame(round(M_SBRC_UKB$prs,3))

write.table(prs_PRSCS_ORG,file="b4u_bmi_prs_PRSCS_ORG.txt",row.names=FALSE,
    col.names=FALSE,quote=FALSE)
write.table(prs_PRSCS_ROB,file="b4u_bmi_prs_PRSCS_ROB.txt",row.names=FALSE,
    col.names=FALSE,quote=FALSE)
write.table(prs_PRSCS_REC,file="b4u_bmi_prs_PRSCS_REC.txt",row.names=FALSE,
    col.names=FALSE,quote=FALSE)
write.table(prs_SBRC_TGP,file="b4u_bmi_prs_SBRC_TGP.txt",row.names=FALSE,
    col.names=FALSE,quote=FALSE)
write.table(prs_SBRC_UKB,file="b4u_bmi_prs_SBRC_UKB.txt",row.names=FALSE,
    col.names=FALSE,quote=FALSE)
```

Now, create an archive with all the produced files to be uploaded:

```bash
SITE="YOUR_PARTNER_NAME_HERE_EG_VUA"

tar -czvf $SITE".tar.gz" b4u_bmi_prs_*.txt
```

Upload the file `SITE.tar.gz` [here](https://drive.google.com/drive/folders/1q3Vl0mfvJt4TXKCJUo4KU7vGiaFo0iCj?usp=drive_link)

Thank you!

**Important note**: There are tiny differences (range -10e-6, 10e-6 for HUA) 
between the two approaches (single file vs chromosome-wise files) when 
calculating PRS with PLINK. These differences are attributed to floating point 
rounds and the order of calculations in PLINK. In any case, these differences
should not have any significant impact in actual PRS calculations or genetic
risk.

