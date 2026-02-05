BETTER4U Task 3.7 - LOO PRS evaluation for BMI _only_
================================================================================

## Authors

**Panagiotis Moulos** <sup>1,2</sup><br/>
<sup>1</sup>Harokopio University of Athens, <sup>2</sup>BSRC Alexander Fleming

### Contact details

* Panagiotis Moulos (moulos@fleming.gr)
* Jon Anders Eriksson (anders.eriksson@ut.ee)

# Introduction

This document explains the PRS evaluation process after the discussed 
Leave-One-Out PRS derivation analysis for BETTER4U cohorts. Each partner should
use the provided PRS score files that have been produced with the respective
partner cohort left-out from the derivation process. In practice, this means 
that partner **X** will use the genotype information of cohort **X cohort** to
produce PLINK score profiles using the provided PRS score files (containing 
SNPs, risk alleles and weights/betas) that have been created with meta-analyzed 
summary statistics calculated _without_ the presence of **X cohort**. 

As of 12/12/2025, LOO summary statistics are available for the cohorts below.
The names indicate that they are _left out_ from the meta-analysis, i.e. BIB
means that meta-analysis was conducted with BIONIC, HMGU, HUA, MUW, TAUH, UH,
UTARTU, VUA:

- BIB
- BIONIC
- HMGU
- HUA
- MUW
- TAUH
- UH (FTC)
- UTARTU
- VUA

The data that will be collected from each left-out cohort can be summarized to:

* Calculation of performance metrics for three PRS versions
  - A GCTB SBayesRC PRS derived with a custom LD panel from 1000 genomes EUR 
    samples
  - A GCTB SBayesRC PRS derived with the built-in UKB LD panel
  - A PRS-CS PRS derived with the built-in EUR 1000 genomes panel and a BIM file 
    constructed from the BETTER4U BMI summary statistics

For the completion of the above tasks, a set of files will be distributed along
with instructions of how to use them. Please make sure that you maintain your
cohort data either as PLINK BED files or Oxford BGEN files, either as a whole or
split per chromosome.

## Required software

* [PLINK 1.90](https://www.cog-genomics.org/plink/)
* [PLINK 2.00](https://www.cog-genomics.org/plink/2.0/) (for BGEN instead of BED input)
* [R](https://www.r-project.org/)
* [gdown](https://github.com/wkentaro/gdown) (optional)
* [This](https://github.com/moulos-lab/better4u/blob/main/03_prs_derivation/evalfuns.R) set of R scripts

## Requested work description

In the following, we assume that the prefix name of your cohort is `COHORT`. For
the `COHORT` we suppose that one of the following applies:

a1. There is one triplet (BED, BIM, FAM) of PLINK files for all the variants in
your cohort named:

```
COHORT.bed
COHORT.bim
COHORT.fam
```

a2. There is one pair (BGEN, sample) of Oxford BGEN files for all the variants 
in your cohort named:

```
COHORT.bgen
COHORT.sample
```

b1. There is one triplet of PLINK files per chromosome, named:

```
COHORT_chrZ.bed
COHORT_chrZ.bim
COHORT_chrZ.fam
```

b2. There is one pair of BGEN files per chromosome, named:

```
COHORT_chrZ.bgen
COHORT_chrZ.sample
```

where `Z = 1..22`. **Note** that the separator between `COHORT` and `chrZ`
should be an underscore (`_`). Altough in the accompanying scripts the separator
can be specified, we suggest to keep it to undescore (`_`). You can explore the
`sanitizePrs()` or `evalPrs()` functions in `evalfuns.R`.

## Calculation of PRS for BMI

### Run setup

Setup the working folder. This can be a folder of your preference in your
system, e.g. `/home/myself/better4u/prs`:

```bash
# The prepending export is important, please do not omit!
export WORKSPACE=/home/myself/better4u/prs
mkdir -p $WORKSPACE && cd $WORKSPACE

git clone https://github.com/moulos-lab/better4u.git
```

### PRS score file retrieval

Download the 3 PRS score files corresponding to your left-out cohort. The 
analyst responsible for each left-out cohort is expected to download a PRS for 
the following cases (3 in total):

1. PRS created with GCTB SBayesRC with custom LD panel from 1000 genomes EUR 
samples
2. PRS created with GCTB SBayesRC with built-in UKB LD panel
3. PRS created with PRS-CS with built-in EUR 1000 genomes panel and a BIM file 
constructed from the BETTER4U BMI summary statistics

The files are named with the following pattern:

- `b4u_bmi_sbrc_tgp_{COHORT}_out.prs` for case (1)
- `b4u_bmi_sbrc_ukb_{COHORT}_out.prs` for case (2)
- `b4u_bmi_prscs_tgp_{COHORT}_out.prs` for case (3)

You can find these files in [this](https://drive.google.com/drive/folders/19unmYal4uqk97jpZKwSSNuaaVyGWSD8n?usp=sharing)
Google Drive link. Each subdirectory is named according to the left-out cohort.
YOu should download and place them inside the working directory. Files where
SNPs are named based on `CHR:POS:REF:ALT` are also available. For your 
convenience, we also provide `gdown` download links, if you have 
[gdown](https://github.com/wkentaro/gdown) installed, otherwise just download
the files and place to the working folder:

#### BIB

```bash
# b4u_bmi_sbrc_tgp_BIB_out.prs
gdown 10F6O-SFC-5FLBC_2CrSCkTsQcrIlH2vI --output b4u_bmi_sbrc_tgp_BIB_out.prs

# b4u_bmi_sbrc_ukb_BIB_out.prs
gdown 1vRF62465gqVR1-XzwYDcVMZQ7B681KgS --output b4u_bmi_sbrc_ukb_BIB_out.prs

# b4u_bmi_prscs_tgp_BIB_out.prs
gdown 11OUijyhx5bh-Fu-yB43zu5NiQsi6ujag --output b4u_bmi_prscs_tgp_BIB_out.prs
```

**or if you prefer `CHR:POS:REF:ALT` SNP ID format:**

```bash
# b4u_bmi_sbrc_tgp_BIB_out.prs
gdown 11H1_kwPiuuBA8y2KinmxxEgfDRYCWRrU --output b4u_bmi_sbrc_tgp_BIB_out.prs

# b4u_bmi_sbrc_ukb_BIB_out.prs
gdown 1w3KX4M0zHyhcITblwh_Amo1L-MTw3x2Y --output b4u_bmi_sbrc_ukb_BIB_out.prs

# b4u_bmi_prscs_tgp_BIB_out.prs
gdown 1g-o6tizFvlqO6R3ZvVzRem0Z9jJ13RG_ --output b4u_bmi_prscs_tgp_BIB_out.prs
```

#### BIONIC

```bash
# b4u_bmi_sbrc_tgp_BIONIC_out.prs
gdown 1ONrJFrHOinBjfIzj4-V2bzZAoC44MByT --output b4u_bmi_sbrc_tgp_BIONIC_out.prs

# b4u_bmi_sbrc_ukb_BIONIC_out.prs
gdown 1RGPJtzscWcK91IluE6DAMpLP7iEvbzOk --output b4u_bmi_sbrc_ukb_BIONIC_out.prs

# b4u_bmi_prscs_tgp_BIONIC_out.prs
gdown 1qCxx1UKafra343-1P6XJEi7HhnMw1el2 --output b4u_bmi_prscs_tgp_BIONIC_out.prs
```

**or if you prefer `CHR:POS:REF:ALT` SNP ID format:**

```bash
# b4u_bmi_sbrc_tgp_BIONIC_out.prs
gdown 1nb29j_7KtU8HIDtz_AC64IT_DLH0M1Rd --output b4u_bmi_sbrc_tgp_BIONIC_out.prs

# b4u_bmi_sbrc_ukb_BIONIC_out.prs
gdown 1M7Ok2iHypMM-z923O_Wcq-FPNnJ-h8rU --output b4u_bmi_sbrc_ukb_BIONIC_out.prs

# b4u_bmi_prscs_tgp_BIONIC_out.prs
gdown 1lY5qPmanLVPR_ybo7B9J8vuAgGU5Slww --output b4u_bmi_prscs_tgp_BIONIC_out.prs
```

#### HMGU

```bash
# b4u_bmi_sbrc_tgp_HMGU_out.prs
gdown 13n6mceDou3p3LQPLi1BGJqi_sL3qkOsl --output b4u_bmi_sbrc_tgp_HMGU_out.prs

# b4u_bmi_sbrc_ukb_HMGU_out.prs
gdown 1GuI2AheG6GRpmEJCSz7Q-kUc71tYkqMq --output b4u_bmi_sbrc_ukb_HMGU_out.prs

# b4u_bmi_prscs_tgp_HMGU_out.prs
gdown 1p4HxWXWNea5GkmbkPf-YstLvX8-6zOvZ --output b4u_bmi_prscs_tgp_HMGU_out.prs
```

**or if you prefer `CHR:POS:REF:ALT` SNP ID format:**

```bash
# b4u_bmi_sbrc_tgp_HMGU_out.prs
gdown 1zzvA2bb-ixDa-z6f3H1dxl4zPO00H0JY --output b4u_bmi_sbrc_tgp_HMGU_out.prs

# b4u_bmi_sbrc_ukb_HMGU_out.prs
gdown 1FjLV1yaiot0e1N7q3x3QD7fMCkJw8cCI --output b4u_bmi_sbrc_ukb_HMGU_out.prs

# b4u_bmi_prscs_tgp_HMGU_out.prs
gdown 12O-mB-VIpAiRFDMRoWO68-BH196-h74Z --output b4u_bmi_prscs_tgp_HMGU_out.prs
```

#### HUA

```bash
# b4u_bmi_sbrc_tgp_HUA_out.prs
gdown 1CLY_hSDPxkT4Mo-51AEmkrKvzl5uCiYn --output b4u_bmi_sbrc_tgp_HUA_out.prs

# b4u_bmi_sbrc_ukb_HUA_out.prs
gdown 15KwAELxNZZMJvtYpCtGm4BTgTiRJn2PE --output b4u_bmi_sbrc_ukb_HUA_out.prs

# b4u_bmi_prscs_tgp_HUA_out.prs
gdown 1rD6-rQHXUhrR4feyN5ktBhpHa908wj7e --output b4u_bmi_prscs_tgp_HUA_out.prs
```

**or if you prefer `CHR:POS:REF:ALT` SNP ID format:**

```bash
# b4u_bmi_sbrc_tgp_HUA_out.prs
gdown 1K4Tmod8IirKr-qqmZOBLlOw2nR51h8R8 --output b4u_bmi_sbrc_tgp_HUA_out.prs

# b4u_bmi_sbrc_ukb_HUA_out.prs
gdown 1y9dEf8ApMemSUnKdeA5FWN2C_y3Uf_7y --output b4u_bmi_sbrc_ukb_HUA_out.prs

# b4u_bmi_prscs_tgp_HUA_out.prs
gdown 1Dsya9f6UfHRhJPPktqhoP7v6w5ko8bDb --output b4u_bmi_prscs_tgp_HUA_out.prs
```

#### MUW

```bash
# b4u_bmi_prscs_tgp_MUW_out.prs
gdown 1glJ1f_yx0hMJN4b3X8Fjt2hI0oR_SYrE --output b4u_bmi_prscs_tgp_MUW_out.prs

# b4u_bmi_sbrc_ukb_MUW_out.prs
gdown 1zCONruHZAZbYOnqjbS4AsecB_7hY16Sh --output b4u_bmi_sbrc_ukb_MUW_out.prs

# b4u_bmi_prscs_tgp_MUW_out.prs
gdown 1Ra7qya5iUvih2EizKBiinVjTdADYF105 --output b4u_bmi_prscs_tgp_MUW_out.prs
```

**or if you prefer `CHR:POS:REF:ALT` SNP ID format:**

```bash
# b4u_bmi_sbrc_tgp_MUW_out.prs
gdown 1ByF1QSjcc0QcMp2JrMzajCfumUcbZeKV --output b4u_bmi_sbrc_tgp_MUW_out.prs

# b4u_bmi_sbrc_ukb_MUW_out.prs
gdown 1C1n0iHC3PQJgi8u-1dIQ3-t5EdeLTj5k --output b4u_bmi_sbrc_ukb_MUW_out.prs

# b4u_bmi_prscs_tgp_MUW_out.prs
gdown 1ezv1U0fA_KVIUzaJdLvvY1uNlE3zwAR2 --output b4u_bmi_prscs_tgp_MUW_out.prs
```

#### TAUH

```bash
# b4u_bmi_sbrc_tgp_TAUH_out.prs
gdown 1DQ_rBwDObkT058oHeRZDIaf_VtItPs7G --output b4u_bmi_sbrc_tgp_TAUH_out.prs

# b4u_bmi_sbrc_ukb_TAUH_out.prs
gdown 1pQy6MeNmYwHaw8MayOb4aL8ed53uzRTd --output b4u_bmi_sbrc_ukb_TAUH_out.prs

# b4u_bmi_prscs_tgp_TAUH_out.prs
gdown 1EwW0n1mCszONjpk5G_n64oAfhXY1naJL --output b4u_bmi_prscs_tgp_TAUH_out.prs
```

**or if you prefer `CHR:POS:REF:ALT` SNP ID format:**

```bash
# b4u_bmi_sbrc_tgp_TAUH_out.prs
gdown 1uHlPBguL3zEgyBUdiXfImwVMqdFfGyiS --output b4u_bmi_sbrc_tgp_TAUH_out.prs

# b4u_bmi_sbrc_ukb_TAUH_out.prs
gdown 18f8hwMaWCtjjgm4IcCchn93VSpk_MIik --output b4u_bmi_sbrc_ukb_TAUH_out.prs

# b4u_bmi_prscs_tgp_TAUH_out.prs
gdown 1AntPLNkh2n-8kz6hxcClCI3thnsZtOkS --output b4u_bmi_prscs_tgp_TAUH_out.prs
```

#### UH

```bash
# b4u_bmi_sbrc_tgp_UH_out.prs
gdown 1B-a2HkgqMuMtQxlXjdXj2Z6VJDopoHjV --output b4u_bmi_sbrc_tgp_UH_out.prs

# b4u_bmi_sbrc_ukb_UH_out.prs
gdown 1gvT2dtRHNnY-igY5V9YL-b3kt2A1FlMD --output b4u_bmi_sbrc_ukb_UH_out.prs

# b4u_bmi_prscs_tgp_UH_out.prs
gdown 1_oRwvX6t_45Ssrby5kMuP3_4_W00O4kz --output b4u_bmi_prscs_tgp_UH_out.prs
```

**or if you prefer `CHR:POS:REF:ALT` SNP ID format:**

```bash
# b4u_bmi_sbrc_tgp_UH_out.prs
gdown 1B_5of5KGN655x6f2h5MujevaoDius8xs --output b4u_bmi_sbrc_tgp_UH_out.prs

# b4u_bmi_sbrc_ukb_UH_out.prs
gdown 1T8-o3cyYavUw43VIJdU64epyPPNvO2N4 --output b4u_bmi_sbrc_ukb_UH_out.prs

# b4u_bmi_prscs_tgp_UH_out.prs
gdown 1O28jpxJk11MCLHk98HXokfK0eqLleZDR --output b4u_bmi_prscs_tgp_UH_out.prs
```

#### UTARTU

```bash
# b4u_bmi_sbrc_tgp_UTARTU_out.prs
gdown 1gVxGrezId1qlRyzOQaUWa6LsAR64rVOn --output b4u_bmi_sbrc_tgp_UTARTU_out.prs

# b4u_bmi_sbrc_ukb_UTARTU_out.prs
gdown 1f63ZSufCM28Xx9YCrfH_PLFmKOXPHfF_ --output b4u_bmi_sbrc_ukb_UTARTU_out.prs

# b4u_bmi_prscs_tgp_UTARTU_out.prs
gdown 1m5MVh04xtdD6-oKNGMgBpH7tNBisJohV --output b4u_bmi_prscs_tgp_UTARTU_out.prs
```

**or if you prefer `CHR:POS:REF:ALT` SNP ID format:**

```bash
# b4u_bmi_sbrc_tgp_UTARTU_out.prs
gdown 1jetJ-nxZgW4GqMcMSrw4J3pB_AeobOUD --output b4u_bmi_sbrc_tgp_UTARTU_out.prs

# b4u_bmi_sbrc_ukb_UTARTU_out.prs
gdown 1n9jfZxIKM5LS1v8KgYvOgv8_xzcB74oh --output b4u_bmi_sbrc_ukb_UTARTU_out.prs

# b4u_bmi_prscs_tgp_UTARTU_out.prs
gdown 1dLQNCwfV5YdStnyTbZWSfHegTn-2XwHI --output b4u_bmi_prscs_tgp_UTARTU_out.prs
```

#### VUA

```bash
# b4u_bmi_sbrc_tgp_VUA_out.prs
gdown 1_qbOj4r8j3SyHJjR36WHEBViNJjBRwex --output b4u_bmi_sbrc_tgp_VUA_out.prs

# b4u_bmi_sbrc_ukb_VUA_out.prs
gdown 1PE935x4AOGCj6zBS-87ZHzXdrNzNO7PD --output b4u_bmi_sbrc_ukb_VUA_out.prs

# b4u_bmi_prscs_tgp_VUA_out.prs
gdown 1noxG7PXk4HAyLJ8jhoW89L2b9NjcdIcG --output b4u_bmi_prscs_tgp_VUA_out.prs
```

**or if you prefer `CHR:POS:REF:ALT` SNP ID format:**

```bash
# b4u_bmi_sbrc_tgp_VUA_out.prs
gdown 10CRMbyNqynTMf2f2Y8MdG9hyJRUZTCNS --output b4u_bmi_sbrc_tgp_VUA_out.prs

# b4u_bmi_sbrc_ukb_VUA_out.prs
gdown 1l3IcA_Y6lfKZ4XYpFuyMPDVMchYPonPN --output b4u_bmi_sbrc_ukb_VUA_out.prs

# b4u_bmi_prscs_tgp_VUA_out.prs
gdown 1KENW5OJ0DqBUM1i7breYEt4PHMg3aO9v --output b4u_bmi_prscs_tgp_VUA_out.prs
```

### Evaluation and PRS calculation

Calculate baseline PRSs based on the downloaded files for each left-out cohort. 
We assume that all the covariate files used in fastGWA analysis are placed in 
the file `covariates.txt` and is present in `$WORKSPACE`. The same covariates as
the ones used with GCTA `--fastGWA` should be used. Note that the trait should a
lso be attached to this `covariates.txt` file. Then, within R:

```R
WORKSPACE <- Sys.getenv("WORKSPACE")

source(file.path(WORKSPACE,"better4u","03_prs_derivation","evalfuns.R"))

# Genetic data file base (whether single or per chromosome, for PLINK or BGEN)
# e.g. COHORT.bed COHORT.bim COHORT.fam
# e.g. COHORT.bgen COHORT.bim COHORT.sample
genoBase <- "COHORT"

# Covariates file (sex, age, age^2, PCs) and trait
covFile <- "covariates.txt"

# Assuming this is the BMI name in covariates, otherwise adjust accordingly
trait <- "bmi"

# Your cohort name, e.g. HUA
cohort <- "YOUR_COHORT_NAME"

# Then, with the downloaded PRS score files

# GCTB SBayesRC with TGP panel
prsFile_SBRC_TGP <- paste0("b4u_bmi_sbrc_tgp_",cohort,"_out.prs")

# GCTB SBayesRC with built-in LD panel
prsFile_SBRC_UKB <- paste0("b4u_bmi_sbrc_ukb_",cohort,"_out.prs")

# PRS-CS with built-in EUR 1000 genomes panel
prsFile_PRSCS_TGP <- paste0("b4u_bmi_prscs_tgp_",cohort,"_out.prs")

# Sanitize the PRS scores (allele flip checking, coverage, etc)
# Single PLINK file:
sanFile_SBRC_TGP <- sanitizePrs(prsFile_SBRC_TGP,genoBase,from="sbayesrc")
sanFile_SBRC_UKB <- sanitizePrs(prsFile_SBRC_UKB,genoBase,from="sbayesrc")
sanFile_PRSCS_TGP <- sanitizePrs(prsFile_PRSCS_TGP,genoBase,from="ready")

# PLINK files per chromosome - if required
#sanFile_SBRC_TGP <- sanitizePrs(prsFile_SBRC_TGP,genoBase,perChr=TRUE,
#   from="sbayesrc",rc=0.2)
#sanFile_SBRC_UKB <- sanitizePrs(prsFile_SBRC_UKB,genoBase,perChr=TRUE,
#   from="sbayesrc",rc=0.2)
#sanFile_PRSCS_TGP <- sanitizePrs(prsFile_PRSCS_TGP,genoBase,perChr=TRUE,
#   from="ready",rc=0.2)

# BGEN file (PLINK 2.0 required in path!):
# sanFile_SBRC_TGP <- sanitizePrs(prsFile_SBRC_TGP,genoBase,bgen=TRUE,
#   from="sbayesrc")
# sanFile_SBRC_UKB <- sanitizePrs(prsFile_SBRC_UKB,genoBase,bgen=TRUE,
#   from="sbayesrc")
# sanFile_PRSCS_TGP <- sanitizePrs(prsFile_PRSCS_TGP,genoBase,bgen=TRUE,
#   from="ready")

# ΒGEN files per chromosome - if required
#sanFile_SBRC_TGP <- sanitizePrs(prsFile_SBRC_TGP,genoBase,perChr=TRUE,
#   bgen=TRUE,from="sbayesrc",rc=0.2)
#sanFile_SBRC_UKB <- sanitizePrs(prsFile_SBRC_UKB,genoBase,perChr=TRUE,
#   bgen=TRUE,from="sbayesrc",rc=0.2)
#sanFile_PRSCS_TGP <- sanitizePrs(prsFile_PRSCS_TGP,genoBase,perChr=TRUE,
#   bgen=TRUE,from="ready",rc=0.2)

# Three files are written:
# b4u_bmi_sbrc_tgp.prs.san
# b4u_bmi_sbrc_ukb.prs.san
# b4u_bmi_prscs_tgp.prs.san

# Then feed to the PRS/metrics function calculation
# Single PLINK file:
M_SBRC_TGP <- evalPrs(sanFile_SBRC_TGP,covFile,trait,genoBase)
M_SBRC_UKB <- evalPrs(sanFile_SBRC_UKB,covFile,trait,genoBase)
M_PRSCS_TGP <- evalPrs(sanFile_PRSCS_TGP,covFile,trait,genoBase)

# Or if you have BGEN + SAMPLE file instead of PLINK files
#M_SBRC_TGP <- evalPrs(sanFile_SBRC_TGP,covFile,trait,genoBase,bgen=TRUE)
#M_SBRC_UKB <- evalPrs(sanFile_SBRC_UKB,covFile,trait,genoBase,bgen=TRUE)
#M_PRSCS_TGP <- evalPrs(sanFile_PRSCS_TGP,covFile,trait,genoBase,bgen=TRUE)

# PLINK files per chromosome - if required
#M_SBRC_TGP <- evalPrs(sanFile_SBRC_TGP,covFile,trait,genoBase,
#   perChr=TRUE,rc=0.2)
#M_SBRC_UKB <- evalPrs(sanFile_SBRC_UKB,covFile,trait,genoBase,
#   perChr=TRUE,rc=0.2)
#M_PRSCS_TGP <- evalPrs(sanFile_PRSCS_TGP,covFile,trait,genoBase,
#   perChr=TRUE,rc=0.2)

# Or BGEN + SAMPLE files per chromosome - if required
#M_SBRC_TGP <- evalPrs(sanFile_SBRC_TGP,covFile,trait,genoBase,
#   bgen=TRUE,perChr=TRUE,rc=0.2)
#M_SBRC_UKB <- evalPrs(sanFile_SBRC_UKB,covFile,trait,genoBase,
#   bgen=TRUE,perChr=TRUE,rc=0.2)
#M_PRSCS_TGP <- evalPrs(sanFile_PRSCS_TGP,covFile,trait,genoBase,
#   bgen=TRUE,perChr=TRUE,rc=0.2)

# Some additional values regarding coverage
N_SBRC_TGP <- as.numeric(countLines(prsFile_SBRC_TGP)) - 1
N_SBRC_UKB <- as.numeric(countLines(prsFile_SBRC_UKB)) - 1
N_PRSCS_TGP <- as.numeric(countLines(prsFile_PRSCS_TGP)) - 1

# Gather metrics to share and PRS (just values, no ids of any kind)
metrics <- data.frame(
    sbayes_tgp=formatC(M_SBRC_TGP$metrics,digits=6),
    sbayes_ukb=formatC(M_SBRC_UKB$metrics,digits=6),
    prscs_tgp=formatC(M_PRSCS_TGP$metrics,digits=6)
)

# Additonal data
add <- as.data.frame(rbind(
    as.integer(c(N_SBRC_TGP,N_SBRC_UKB,N_PRSCS_TGP)),
    round(100*as.numeric(metrics[nrow(metrics),])/
        c(N_SBRC_TGP,N_SBRC_UKB,N_PRSCS_TGP),digits=2)
))
rownames(add) <- c("snps_total","coverage")
colnames(add) <- names(metrics)

# Final metrics data frame
finalMetrics <- rbind(metrics,add)

outFile <- paste0("b4u_bmi_prs_metrics_",cohort,".txt")
write.table(finalMetrics,file=outFile,sep="\t",quote=FALSE,col.names=NA)
```

Upload the file `b4u_bmi_prs_metrics_{YOUR_COHORT_NAME}.txt` [here](https://drive.google.com/drive/folders/1ifphpUZJ3aM_mRq29DfHr6fFapClj7Ls?usp=sharing)

### PRS values for federated learning

The actual PRS values must be calculated and included in the federated learning
of WP5. To this end, export **anonymized** and **truncated** PRS values in a
single column text file as follows:

```R
# Write the actual PRS values for each version
prs_SBRC_TGP <- as.data.frame(round(M_SBRC_TGP$prs,3))
prs_SBRC_UKB <- as.data.frame(round(M_SBRC_UKB$prs,3))
prs_PRSCS_TGP <- as.data.frame(round(M_PRSCS_TGP$prs,3))

write.table(prs_SBRC_TGP,file="b4u_bmi_prs_SBRC_TGP.txt",row.names=FALSE,
    col.names=FALSE,quote=FALSE)
write.table(prs_SBRC_UKB,file="b4u_bmi_prs_SBRC_UKB.txt",row.names=FALSE,
    col.names=FALSE,quote=FALSE)
write.table(prs_PRSCS_TGP,file="b4u_bmi_prs_PRSCS_ORG.txt",row.names=FALSE,
    col.names=FALSE,quote=FALSE)
```

The file `b4u_bmi_prs_SBRC_UKB.txt` should be used for the federated learning
run.

**Thank you in advance**

### Important note

There are tiny differences (range -10e-6, 10e-6 for HUA) between the two 
approaches (single file vs chromosome-wise files) when calculating PRS with 
PLINK. These differences are attributed to floating point rounds and the order 
of calculations in PLINK. In any case, these differences should not have any
significant impact in actual PRS calculations or genetic risk.
