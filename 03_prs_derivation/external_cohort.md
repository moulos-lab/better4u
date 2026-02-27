Evaluation of BETTER4U PGS for Weight Change in External Cohorts
================================================================================

## Authors

**Panagiotis Moulos** <sup>1,2</sup><br/>
<sup>1</sup>Harokopio University of Athens, <sup>2</sup>BSRC Alexander Fleming

### Contact details

* Panagiotis Moulos (moulos@fleming.gr)
* Jon Anders Eriksson (anders.eriksson@ut.ee)
* René Pool (r.pool@vu.nl)
* Dorret Boomsma (di.boomsma@vu.nl)

# Introduction

The [BETTER4U](https://better4u.eu/) consortium recently produced a new GWAS for 
rate of weight change that has identified a large number of novel weight change 
loci (manuscript in preparation). The aim of the present document is to provide
instructions and scripts regarding:

- Preparation of the phenotypic file where the validation model will operate
- Calculation of the weight change phenotype with the same methodology as in 
BETTER4U GWAS
- Calculation of the PGS given BETTER4U PGS files
- Caclulation of PGS evaluation metrics

The population of interest concerns adults. The rest of the document describes
in detail the various steps that have to be followed.

# Implementation

## Required software

* [PLINK 1.90](https://www.cog-genomics.org/plink/)
* [PLINK 2.00](https://www.cog-genomics.org/plink/2.0/)
* [KING](https://www.kingrelatedness.com/Download.shtml)
* [FlashPCA2](https://github.com/gabraham/flashpca)
* [R](https://www.r-project.org/)
* [gdown](https://github.com/wkentaro/gdown) (optional)
* [This](https://github.com/moulos-lab/better4u/blob/main/03_prs_derivation/evalfuns.R) set of R scripts

## Preparation of phenotypic file

*Template scripts that perform the transformations described below are provided with this document*

A simple text tab-delimited file (`pheno.txt`) should be prepared with the 
following phenotypic information:

| #column | phenotype | description |
| --- | --- | --- |
| 1 | FID | Family ID |
| 2 | IID | Individual ID |
| 3 | Sex | Biological sex 1 (male), 2 (female) encoded |
| 4 | Age | Age at measurement (if relevant, with decimals to distinguish measurements within the same year) |
| 5 | Height | Height at measurement (in units of m) |
| 6 | Weight | Body mass at measurement (in units of kg) |
| 7 | BMI | BMI at measurement (in units of kg/m<sup>2</sup>) |

For individuals with measurements at multiple time points, please combine the 
values for Age, Height, Weight and BMI into single strings with values from each
measurement separated by semicolon (`;`). For example, if individual S001 has 
three measurements:

| #Measurement | Sex | Age | Height | Weight | BMI |
| --- | --- | --- | --- | --- | --- |
| 1 | 0 | 24 | 1.93 | 113.5 | 30.47 |
| 2 | 0 | 28 | 1.93 | 111.5 | 29.92 |
| 3 | 0 | 32 | 1.93 | 109.9 | 29.5 |

Then the file will be formatted as:

| FID | IID | Sex | Age | Height | Weight | BMI |
| --- | --- | --- | --- | --- | --- | --- |
| S001 | S001 | 0 | 24;28;32 | 1.93;1.93;1.93 | 113.5;111.5;109.9 | 30.47;29.92;29.50 |

**Important considerations** 

- Please exclude measurements occurring less than 12 months after the previous 
measurement.
- The order of measurements must be the same in all columns (but do not 
necessarily have to be arranged in order of increasing time).
- Missing values: if not all variables are measured at a given time point, 
please leave the entry blank, or put the string NA for missing value. In the 
example above, if height was only measured at the first time point, we would 
have the line

| FID | IID | Sex | Age | Height | Weight | BMI |
| --- | --- | --- | --- | --- | --- | --- |
| S001 | S001 | 0 | 24;28;32 | 1.93;; | 113.5;111.5;109.9 | 30.47;29.92;29.50 |

The number of semicolons must always be one less than the number of data points 
for the individual. Alternatively, using NA instead of a blank entry:

| FID | IID | Sex | Age | Height | Weight | BMI |
| --- | --- | --- | --- | --- | --- | --- |
| S001 | S001 | 0 | 24;28;32 | 1.93;NA;NA | 113.5;111.5;109.9 | 30.47;29.92;29.50 |

## Calculation of summary statistics for BMI and Weight change

This section describes the process of calculating BMI and weight change slopes
and provides related scripts. The resulting file will have the following columns
and format (tab-delimited):

| #column | name | description |
| --- | --- | --- |
| 1 | FID | Family ID |
| 2 | IID | Individual ID |
| 3 | N | Number of measurements (N >= 2) |
| 4 | M_Age | Mean age at measurement |
| 5 | S_Age | Standard error of age at measurement |
| 6 | M_BMI | Mean BMI |
| 7 | S_BMI | Standard error of BMI |
| 8 | BMI_alpha | Intercept from linear regression `BMI_t ~ alpha + beta * (Age_t - Age_start)` where `Age_start` is the lower bound on age for the category |
| 9 | BMI_beta | Slope of the linear regression for BMI |
| 10 | M_Weight | Mean weight |
| 11 | S_Weight | Standard error of weight |
| 12 | Weight_alpha | Intercept from linear regression `Weight_t ~ alpha + beta * (Age_t - Age_start)` where `Age_start` is the lower bound on age for the category |
| 13 | Weight_beta | Slope of the linear regression for weight |
| 14 | Height_mean | Mean value of reported height measurements |

The following R script will read the file `pheno.txt` prepared above and 
calculate the aforementioned summary statistics:

```r
# Read in the prepared file
d  <- read.delim("pheno.txt")
nr <- nrow(d)

# Init output variables
d$N <- d$M_Age <- d$S_Age <- d$M_BMI <- d$S_BMI <- d$BMI_Alpha <- 
    d$BMI_Beta <- d$M_Weight <- d$S_Weight <- d$Weight_Alpha <- 
    d$Weight_Beta <- d$Height_mean <- rep(0,nr)

for (k in 1:nr) {
    # Split Age, BMI and Weight strings into vectors
    age_t <- suppressWarnings(as.numeric(strsplit(d$Age[k],";")[[1]]))
    bmi <- suppressWarnings(as.numeric(strsplit(d$BMI[k],";")[[1]]))
    weight <- suppressWarnings(as.numeric(strsplit(d$Weight[k],";")[[1]]))
    
    # Calculate mean height from available values
    if (is.numeric(d$Height[k])) # One value only
        d$Height_mean[k] <- d$Height[k]
    else {
        h <- suppressWarnings(as.numeric(strsplit(d$Height[k],";")[[1]]))
        d$Height_mean[k] <- mean(h,na.rm=TRUE)
    }

    # Select time points with BMI information and age in the right interval
    use <- !is.na(age_t) & !is.na(bmi) & !is.na(weight) & age_t > 18 #& age_t <= 65

    # Subest data to eligible time points and offset so that Alpha is estimate 
    # of BMI/Weight at age 18
    age_t <- age_t[use] - 18
    bmi <- bmi[use]
    weight <- weight[use]

    d$N[k] <- length(age_t)
    if (d$N[k] >= 2) {
        d$M_Age[k] <- mean(age_t);
        d$M_BMI[k] <- mean(bmi);
        d$M_Weight[k] <- mean(weight);
        
        # Unbiased estimator of standard error of age in the sample
        d$S_Age[k] <- sqrt(sum((age_t-mean(age_t))^2)/(d$N[k]-1));
        # Unbiased estimator of standard error of BMI in the sample
        d$S_BMI[k] <- sqrt(sum((bmi-mean(bmi))^2)/(d$N[k]-1));
        # Unbiased estimator of standard error of Weight in the sample
        d$S_Weight[k] = sqrt(sum((weight-mean(weight))^2)/(d$N[k]-1)); 
        
        # Models
        f1 <- summary(lm(bmi ~ age_t))
        f2 = summary(lm(weight ~ age_t))
        
        # Estimator of intercept at age 18 - BMI
        d$BMI_Alpha[k] <- f1$coefficients[1,1]
        # Estimator of slope of BMI
        d$BMI_Beta[k] <- f1$coefficients[2,1]
        # Estimator of intercept at age 18 - Weight
        d$Weight_Alpha[k] <- f2$coefficients[1,1]
        # Estimator of slope of BMI
        d$Weight_Beta[k] <- f2$coefficients[2,1]
    }
}

# Save data to file for individuals with >= 2 measurement points
cols_to_keep <- c("FID","IID","Sex","N","Height_mean","M_Age","S_Age","M_BMI",
    "S_BMI","M_Weight","S_Weight","BMI_Alpha","BMI_Beta","Weight_Alpha",
    "Weight_Beta")
d_keep <- d[d$N>=2,cols_to_keep]
write.table(d_keep,"bmi_weight_change_summary_adults.txt",
    sep="\t",quote=FALSE,row.names=FALSE)
```

The exported file `bmi_weight_chage_summary_adults.txt` is the basis for the
steps to follow, that is preparation of genotypic data, principal components and
phenotypic file for regression.

## Preparation of genotypic data with BMI and Weight Change samples

The instructions below are general steps required to derive Principal Components
to include in the PGS validation regression model. Your cohort may already be
processed according to several of the steps below. Please read the following and
perform only the steps you think are necessary given the status of your cohort.

### General assumptions

- The main file format is PLINK 1.90 binary format
- Variants are annotated with dbSNP ids
- Work will be performed with genetic information in autosomes only
- The cohort basename is `COHORT`
- If genotype files are split per chromosome, the naming convention is
`COHORT_chrZ`, where `Z` in `1..22`
- Work is performed in a directory where genetic data files and 
`bmi_weight_chage_summary_adults.txt` are present
- `plink` (and `plink2` if required) are in the system's `$PATH`

In the following, we are working with PLINK 1.90 binary file formats (`.bim + 
.bed + .fam`). Conversion steps between other indicative formats are provided.
To avoid excessive file handling, you can create symlinks to actual data files
and name them `COHORT.*`, for example:

```bash
ln -s /path/to/data/MY_DATA.bed COHORT.bed
ln -s /path/to/data/MY_DATA.bim COHORT.bim
ln -s /path/to/data/MY_DATA.fam COHORT.fam
```

**Note**: The following processes aim at i) identifying individuals with both
phenotypic and genotypic data present, ii) performing genotypic data QC to
derive PCs to include in the regression PGS validation model. QCed genotypic
data may or may not be used for the PGS validation. We suggest using all data
to maximize coverage with the BETTER4U PGS.

### Preparation of individuals whith both genetic and phenotypic data

The file `bmi_weight_change_summary_adults.txt` produced before contains FIDs 
and IIDs from the participants. We will intersect these with the individuals 
present in `COHORT.fam`, or `COHORT.sample` if you have `.bgen` files, or 
`COHORT.psam` if you have PLINK2 `.pgen` files. If files per chromosome are 
present, then we have `COHORT_chrZ.fam` or `COHORT_chrZ.sample` etc.

The following R script will produce a set of individuals to keep.

```r
Rscript \
  -e '{
    # Define basename
    base <- "COHORT"
    # or base <- "COHORT_chr1" if split per chromosome

    # Read-in generated summary statistics
    indiv <- read.delim("bmi_weight_change_summary_adults.txt")

    # Read in samples with genetic data available
    # If your FAM file is tab delimited, use the commented line below
    if (file.exists(paste0(base,".fam")))
        fam <- read.table(paste0(base,".fam"))
        #fam <- read.delim(paste0(base,".fam"),header=FALSE)
    else if (file.exists(paste0(base,".sample"))) {
        fam <- read.table(paste0(base,".sample"),header=TRUE)
        fam <- fam[-1,,drop=FALSE]
    }
    else if (file.exists(paste0(base,".psam")))
        fam <- read.table(paste0(base,".psam"))

    # Assign rownames for easier subsetting
    rownames(fam) <- fam[,2]

    # Intersect and write
    common <- intersect(rownames(fam),indiv$IID)
    write.table(fam[common,1:2,drop=FALSE],file="indiv_keep.txt",
        col.names=FALSE,row.names=FALSE,quote=FALSE)
  }'
```

### Subsetting, file conversion and merging

In this section, we perform file conversions to PLINK 1.90 (if necessary) as 
well as merging multiple chromosomes to a single fileset for downstream 
analysis.

#### Files in PLINK 1.90 format

In this case, only subsetting is necessary. If files are split per chromosome:

```bash
for CHR in `seq 1 22`
do
  plink2 \
    --bfile COHORT_chr${CHR} \
    --keep indiv_keep.txt \
    --snps-only \
    --max-alleles 2 \
    --rm-dup force-first \
    --make-bed \
    --memory 22000 \
    --threads 3 \
    --out COHORT_dbmi_chr${CHR} &
done
```

and then merge:

```bash
for CHR in `seq 1 22`
do
  echo COHORT_dbmi_chr${CHR} >> merge.list
done

plink \
  --merge-list merge.list \
  --make-bed \
  --out COHORT_dbmi
```

If there is only one fileset:

```bash
plink2 \
  --bfile COHORT \
  --keep indiv_keep.txt \
  --snps-only \
  --max-alleles 2 \
  --rm-dup force-first \
  --make-bed \
  --out COHORT_dbmi
```

#### Files in PLINK 2.00 format

Conversion and subsetting need to be performed. If files are split per 
chromosome:

```bash
for CHR in `seq 1 22`
do
  plink2 \
    --pfile COHORT_chr${CHR} \
    --keep indiv_keep.txt \
    --snps-only \
    --max-alleles 2 \
    --rm-dup force-first \
    --make-bed \
    --memory 22000 \
    --threads 3 \
    --out COHORT_dbmi_chr${CHR} &
done
```

and then merge:

```bash
for CHR in `seq 1 22`
do
  echo COHORT_dbmi_chr${CHR} >> merge.list
done

plink \
  --merge-list merge.list \
  --make-bed \
  --out COHORT_dbmi
```

If there is only one fileset:

```bash
plink2 \
  --pfile COHORT \
  --keep indiv_keep.txt \
  --snps-only \
  --max-alleles 2 \
  --rm-dup force-first \
  --make-bed \
  --out COHORT_dbmi
```

#### Files in Oxford BGEN format

Conversion and subsetting need to be performed. If files are split per 
chromosome:

```bash
for CHR in `seq 1 22`
do
  plink2 \
    --bgen COHORT_chr${CHR}.bgen ref-last \
    --sample COHORT_chr${CHR}.sample \
    --keep indiv_keep.txt \
    --snps-only \
    --max-alleles 2 \
    --rm-dup force-first \
    --make-bed \
    --memory 22000 \
    --threads 3 \
    --out COHORT_dbmi_chr${CHR} &
done
```

and then merge:

```bash
for CHR in `seq 1 22`
do
  echo COHORT_dbmi_chr${CHR} >> merge.list
done

plink \
  --merge-list merge.list \
  --make-bed \
  --out COHORT_dbmi
```

If there is only one fileset:

```bash
plink2 \
  --bgen COHORT.bgen ref-last \
  --sample COHORT.sample \
  --keep indiv_keep.txt \
  --snps-only \
  --max-alleles 2 \
  --rm-dup force-first \
  --make-bed \
  --out COHORT_dbmi
```

### Basic SNP and sample filtering

Next we perform filtering with PLINK 1.9 to remove rare variants as well as 
lower sample and call rates. The following variant filters are recommended but
they can be changed according to dataset particularities and known behaviour:

1. Variant call rate: >98%
2. Minor Allele Frequency: >0.05
3. Hardy-Weinberg equilibrium p-value: >10<sup>-9</sup>

The following sample filters are recommended:

1. Sample call rate: > 95%
2. Heterozygosity: median(heterozygosity) &plusmn; 3 &times; IQR
3. KING relatedness: 2nd degree

Below, we sequentially apply variant and sample filters according to widely
accepted [best practices](https://onlinelibrary.wiley.com/doi/10.1002/sim.6605).

First, we calculate heterozygosities:

```bash
plink \
  --bfile COHORT_dbmi \
  --out COHORT_dbmi --het
```

Then, we create a file with sample names with heterozygosities within the limits
of our filter (sample filter #2 above):

```r
Rscript \
  -e '{
    hetData <- read.table("COHORT_dbmi.het",header=TRUE,
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
  --bfile COHORT_dbmi \
  --out COHORT_dbmi_tmp \
  --make-bed \
  --geno 0.02 \
  --maf 0.05 \
  --hwe 0.000000001 \
  --mind 0.05
```

We now create files with variants and samples to *keep* (samples to keep are 
merged with those passing heterozygosity filters):

```r
cut -d" " -f1-2 COHORT_dbmi_tmp.fam > generic_samples_pass.txt
cut -f2 COHORT_dbmi_tmp.bim > generic_variants_pass.txt

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
  --bfile COHORT_dbmi \
  --out COHORT_dbmi_filtered \
  --extract generic_variants_pass.txt \
  --keep all_samples_pass.txt \
  --make-bed
```

Kinship analysis is performed with KING (no LD pruning):

```bash
king \
  -b COHORT_dbmi_filtered.bed \
  --unrelated \
  --degree 2 \
  --prefix individuals_
```

The aforementioned command will create the file `individuals_unrelated.txt` 
which we use with PLINK to create the dataset with related individuals removed
(optional).

```bash
plink \
  --bfile COHORT_dbmi_filtered \
  --keep individuals_unrelated.txt \
  --out COHORT_dbmi_final \
  --make-bed
```

### LD pruning

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
    --bfile COHORT_dbmi_final \
    --indep-pairwise 500 50 $R2 \
    --out ./pruned/COHORT_dbmi_final
    
  TOTAL=$(wc -l ./pruned/*.in | awk '{print $1}')
  
  echo "Parameters --indep-pairwise 500 50 $R2 yield $TOTAL SNPs" >> \
    pruning.log
done
```

The final value of `$R2` can be found by:

```bash
tail -1 pruning.log | awk '{print $(NF-3)}' | sed 's/^\./0./'
```

We use this value to perform the final pruning:

```bash
rm -r ./pruned/*

FINAL_R2=`tail -1 pruning.log | awk '{print $(NF-3)}' | sed 's/^\./0./'`
plink \
  --bfile COHORT_dbmi_final \
  --indep-pairwise 500 50 $FINAL_R2 \
  --out COHORT_dbmi_final_pruned

plink \
  --bfile COHORT_dbmi_final \
  --extract COHORT_dbmi_final_pruned.prune.in \
  --make-bed \
  --out COHORT_dbmi_final_pruned

rm *.nosex *.log
```

### Principal Component Analysis with FlashPca2

Perform PCA using FlashPCA2:

```
flashpca \
  --bfile COHORT_dbmi_final_pruned \
  --ndim 20 \
  --outpc pcs_COHORT_dbmi.txt \
  --outvec eig_COHORT_dbmi.txt \
  --outload loads_COHORT_dbmi.txt
```

## Phenotypic files with Principal Components

Prior to evaluating the BETTER4U Polygenic Scores, we need to construct proper
phenotypic files for regression. These include:

- IID
- Sex
- Age and Age<sup>2</sup>
- BMI change slopes (`BMI_Beta`) for BMI change
- Weight change slopes (`Weight_Beta`) for weight change
- The first 10 PCs

The response variables will be the betas, based on BETTER4U-developed 
methodology.

### Phenotypic file for BMI change

```r
# Read and adjust covariates
covars <- read.delim("bmi_weight_change_summary_adults.txt")
rownames(covars) <- covars$IID
cols_keep <- c("IID","Sex","M_Age","BMI_Beta")
covars <- covars[,cols_keep]
names(covars)[3] <- "Age"
covars$AgeSq <- covars$Age^2

# Attach PCs
pcs <- read.delim("pcs_COHORT_dbmi.txt")
rownames(pcs) <- pcs$IID
covars <- covars[rownames(pcs),,drop=FALSE]
covars10 <- cbind(covars,pcs[,3:12])
write.table(covars10,file="covariates_dbmi_10pcs.txt",sep="\t",quote=FALSE,
    row.names=FALSE)
```

### Phenotypic file for weight change

```r
# Read and adjust covariates
covars <- read.delim("bmi_weight_change_summary_adults.txt")
rownames(covars) <- covars$IID
cols_keep <- c("IID","Sex","M_Age","Weight_Beta")
covars <- covars[,cols_keep]
names(covars)[3] <- "Age"
covars$AgeSq <- covars$Age^2

# Attach PCs
pcs <- read.delim("pcs_COHORT_dbmi.txt")
rownames(pcs) <- pcs$IID
covars <- covars[rownames(pcs),,drop=FALSE]
covars10 <- cbind(covars,pcs[,3:12])
write.table(covars10,file="covariates_dwc_10pcs.txt",sep="\t",quote=FALSE,
    row.names=FALSE)
```

## Assessment of BETTER4U Polygenic Score for Weight Change

The Polygenic Scores were developed using weight change GWAS. However, we 
evaluate both against BMI and weight change. There are two PGS files with
weighted betas derived using [SBayesRC](https://gctbhub.cloud.edu.au/software/gctb/#SBayesRCTutorial)
and [PRS-CS](https://github.com/getian107/PRScs). They can be retrieved from the
following Google Drive links:

1. [SBayesRC PGS](https://drive.google.com/file/d/1hFGxhncwmY4lTMAOtS7wa5Hx7cOyngsC/view?usp=drive_link)
2) [PRS-CS PGS](https://drive.google.com/file/d/15136PvuqV9hpKV-Bc5LZTGv_FUv_Ew1H/view?usp=drive_link)

Alternatively, with `gdown`:

```bash
# PRS-CS PGS
gdown 15136PvuqV9hpKV-Bc5LZTGv_FUv_Ew1H --output b4u_wc_prscs_tgp.prs

# SBayesRC PGS
gdown 1hFGxhncwmY4lTMAOtS7wa5Hx7cOyngsC --output b4u_wc_sbrc_ukb.prs
```

Below we produce evaluation metrics for `COHORT` after basic phenotypic QC. To 
achieve maximum PRS coverage we skip QC in genotypic data.


### Evaluation for BMI change

```r
source("https://raw.githubusercontent.com/moulos-lab/better4u/refs/heads/main/03_prs_derivation/evalfuns.R")

trait <- "BMI_Beta"
genoBase <- "COHORT_dbmi"
covFile <- "covariates_dbmi_10pcs.txt"

# PRS-CS with 1000 genomes LD - original version
prsFile_PRSCS_ORG <- "b4u_wc_prscs_tgp.prs"

## PGS file sanitization (allele flip checks etc.)
# If you have COHORT files per chromosome, use the command below
sanFile_PRSCS_ORG <- sanitizePrs(prsFile_PRSCS_ORG,genoBase,perChr=TRUE,
    from="ready",rc=0.4)
# Otherwise, use the command below
sanFile_PRSCS_ORG <- sanitizePrs(prsFile_PRSCS_ORG,genoBase,from="ready")

## Sanitized PGS evaluation
# If you have COHORT files per chromosome, use the command below
M_PRSCS_ORG <- evalPrs(sanFile_PRSCS_ORG,covFile,trait,genoBase,iidCol=1,
    perChr=TRUE,rc=0.4)
# Otherwise, use the command below
M_PRSCS_ORG <- evalPrs(sanFile_PRSCS_ORG,covFile,trait,genoBase,iidCol=1)

# SBayesRC with built-in UKB LD
prsFile_SBRC_UKB <- "b4u_wc_sbrc_ukb.prs"

## PGS file sanitization (allele flip checks etc.)
# If you have COHORT files per chromosome, use the command below
sanFile_SBRC_UKB <- sanitizePrs(prsFile_SBRC_UKB,genoBase,perChr=TRUE,
    from="ready",rc=0.4)
# Otherwise, use the command below
sanFile_SBRC_UKB <- sanitizePrs(prsFile_SBRC_UKB,genoBase,from="ready")

## Sanitized PGS evaluation
# If you have COHORT files per chromosome, use the command below
M_SBRC_UKB <- evalPrs(sanFile_SBRC_UKB,covFile,trait,genoBase,iidCol=1,
    perChr=TRUE,remFile="removed.txt",rc=0.4)
# Otherwise, use the command below
M_SBRC_UKB <- evalPrs(sanFile_SBRC_UKB,covFile,trait,genoBase,iidCol=1)
    
save(M_PRSCS_ORG,M_SBRC_UKB,file="metrics_dbmi_10pcs.rda")

# Some additional values regarding coverage
N_PRSCS_ORG <- as.numeric(countLines(prsFile_PRSCS_ORG)) - 1
N_SBRC_UKB <- as.numeric(countLines(prsFile_SBRC_UKB)) - 1

# Gather metrics to share and PRS (just values, no ids of any kind)
metrics <- data.frame(
    prscs_tgp=formatC(M_PRSCS_ORG$metrics,digits=6),
    sbayes_ukb=formatC(M_SBRC_UKB$metrics,digits=6)
)

# Additonal data
add <- as.data.frame(rbind(
    as.integer(c(N_PRSCS_ORG,N_SBRC_UKB)),
    round(100*as.numeric(metrics[nrow(metrics),])/
        c(N_PRSCS_ORG,N_SBRC_UKB),digits=2)
))
rownames(add) <- c("snps_total","coverage")
colnames(add) <- names(metrics)

# Final metrics data frame
finalMetrics <- rbind(metrics,add)

write.table(finalMetrics,file="dbmi_prs_metrics_10pcs.txt",sep="\t",
    quote=FALSE,col.names=NA)
```

### Evaluation for Weight change

```r
source("https://raw.githubusercontent.com/moulos-lab/better4u/refs/heads/main/03_prs_derivation/evalfuns.R")

trait <- "Weight_Beta"
genoBase <- "COHORT_dbmi"
covFile <- "covariates_dwc_10pcs.txt"

# PRS-CS with 1000 genomes LD - original version
prsFile_PRSCS_ORG <- "b4u_wc_prscs_tgp.prs"

## PGS file sanitization (allele flip checks etc.)
# If you have COHORT files per chromosome, use the command below
sanFile_PRSCS_ORG <- sanitizePrs(prsFile_PRSCS_ORG,genoBase,perChr=TRUE,
    from="ready",rc=0.4)
# Otherwise, use the command below
sanFile_PRSCS_ORG <- sanitizePrs(prsFile_PRSCS_ORG,genoBase,from="ready")

## Sanitized PGS evaluation
# If you have COHORT files per chromosome, use the command below
M_PRSCS_ORG <- evalPrs(sanFile_PRSCS_ORG,covFile,trait,genoBase,iidCol=1,
    perChr=TRUE,rc=0.4)
# Otherwise, use the command below
M_PRSCS_ORG <- evalPrs(sanFile_PRSCS_ORG,covFile,trait,genoBase,iidCol=1)

# SBayesRC with built-in UKB LD
prsFile_SBRC_UKB <- "b4u_wc_sbrc_ukb.prs"

## PGS file sanitization (allele flip checks etc.)
# If you have COHORT files per chromosome, use the command below
sanFile_SBRC_UKB <- sanitizePrs(prsFile_SBRC_UKB,genoBase,perChr=TRUE,
    from="ready",rc=0.4)
# Otherwise, use the command below
sanFile_SBRC_UKB <- sanitizePrs(prsFile_SBRC_UKB,genoBase,from="ready")

## Sanitized PGS evaluation
# If you have COHORT files per chromosome, use the command below
M_SBRC_UKB <- evalPrs(sanFile_SBRC_UKB,covFile,trait,genoBase,iidCol=1,
    perChr=TRUE,remFile="removed.txt",rc=0.4)
# Otherwise, use the command below
M_SBRC_UKB <- evalPrs(sanFile_SBRC_UKB,covFile,trait,genoBase,iidCol=1)
    
save(M_PRSCS_ORG,M_SBRC_UKB,file="metrics_dwc_10pcs.rda")

# Some additional values regarding coverage
N_PRSCS_ORG <- as.numeric(countLines(prsFile_PRSCS_ORG)) - 1
N_SBRC_UKB <- as.numeric(countLines(prsFile_SBRC_UKB)) - 1

# Gather metrics to share and PRS (just values, no ids of any kind)
metrics <- data.frame(
    prscs_tgp=formatC(M_PRSCS_ORG$metrics,digits=6),
    sbayes_ukb=formatC(M_SBRC_UKB$metrics,digits=6)
)

# Additonal data
add <- as.data.frame(rbind(
    as.integer(c(N_PRSCS_ORG,N_SBRC_UKB)),
    round(100*as.numeric(metrics[nrow(metrics),])/
        c(N_PRSCS_ORG,N_SBRC_UKB),digits=2)
))
rownames(add) <- c("snps_total","coverage")
colnames(add) <- names(metrics)

# Final metrics data frame
finalMetrics <- rbind(metrics,add)

write.table(finalMetrics,file="dwc_prs_metrics_10pcs.txt",sep="\t",
    quote=FALSE,col.names=NA)
```

Please upload the files `metrics_dbmi_10pcs.rda`, `metrics_dwc_10pcs.rda`, `dbmi_prs_metrics_10pcs.txt` and `dwc_prs_metrics_10pcs.txt` in
[this](https://drive.google.com/drive/folders/1OmS3tbUJSNB3OWCU63NLlWB0zYeOL6XK?usp=drive_link)
Google Drive folder and notify the contacts mentioned in the beginning of this
file. If you are not allowed to share the actual polygenic score values 
calculated within the `evalPrs()` function in the code above, please add the
argument `outScores=FALSE`, i.e. `evalPrs(...,outScores=FALSE)`.

Thank you in advance for the cooperation. 
