First Analysis Plan for Genetic Data - WP3
================================================================================

# Introduction

For the start of WP3 of BETTER4U, we will perform genome-wide association 
analyses for body weight and Body Mass Index (BMI). This document focuses on 
tasks **3.2** (*Database of harmonized GWAS summary statistics from consortium
cohorts and publicly available sources*) and **3.4** (*Database of genetic
variants affecting variation in BMI and associated traits, including a list of
genetic loci with pleiotropic effects*). Each individual study will perform data
quality control (QC) and analysis and provide summary results for meta-analysis.
Results files will be deposited to a central repository where additional QC and
the meta-analysis will be performed.

If you have any questions, please email all the following individuals: 

* Jon Anders Eriksson (anders.eriksson@ut.ee)
* Nana Kalafati (nkalafati@gmail.com)
* Panagiotis Moulos (moulos@fleming.gr)
* ...Others


**TIMELINE FOR COMPLETION OF COHORT-SPECIFIC ANALYSES** 

Please submit first results by **MM/DD/2024**

Before starting the analysis, please make sure that:

1. Imputation of genotypes is performed using the HRC r1.1. 2016 reference 
panel. In case your data have not been imputed to this panel yet, please 
proceed as soon as possible.
2. All genotypes are on the correct human genome version build (GRCh37/hg19) and 
on the correct strand (forward). If having used a public imputation server for
HRC r.1.1. 2016, this should hold true.  
3.  Files with genotypic data (arrays/imputed/WGS) are in [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) format.
4. Phenotypic variables are coded based on the [BETTER4U codebook](#).

## Overview 

There are six major components to the WP3 analysis plan: 

1. In this analysis plan, genetic data along with a limited set of variables 
from each cohort will be used to conduct association analyses, namely age, sex, 
body weight and BMI. 
2. Analyses will include both sexes and all age groups. However, analyses will 
be performed separately for individuals over and under 18yrs old, as well as for
males and females.
3. Analyses will not include sex chromosomes (chrX, chrY) or mitochondrial DNA 
(chrM/MT).
4. In studies with more than one measurement points, we ask that you perform the
cross-sectional analyses for each one of your cohort’s measurement points.
5. Template scripts in R will be provided within the [GitHub repository](#) 
linked to project overview XXX. Sharing your work here is highly appreciated. 
These environments are private, so please share your github ID with XXX (XXX).
6. Meta-analysis of the results will be led by the UTARTU team using the 
[MR-MEGA](https://genomics.ut.ee/en/tools#:~:text=Instructions-,MR%2DMEGA,-Introduction) software and methodology. 

## Considerations

Due to the cohorts' diversity in measurements, we would like to ask for the 
analysts to only include data measured using the following measurement units:
1. For BMI, please use kg/m<sup>2</sup> and two decimals (e.g. `23.87`). For BMI
in children and adolescents, please replace BMI with zBMI. You will find 
attached the relevant syntax for calculating zBMI based on 
[Cole and Lobstein et al., 2012](https://pubmed.ncbi.nlm.nih.gov/22715120/) in 
the XXX script. 
2. For weight, please use as unit kg with two decimals. Only perform analyses in
the sample of age greater than 18 years old.
3. For height, please use as unit m in two decimals. Only perform analyses in 
the sample of age greater than 18 years.
4. For sex, please use the code: 1=man, 2=woman.


We thank you for your cooperation.

# Requirements

## General assumptions

1. Quality control steps for samples is going to be performed genome-wide.
2. Quality control steps for variants is going to be performed chromosome-wise
or genome-wide. Certain steps can be performed only genome-wide.
3. Sample/individual names are named reasonably
    1. No special characters, stars, parentheses, dots etc
    2. No language specific characters (e.g. Greek, Swedish, Slavic etc.)
    3. No spaces
4. The main component of the population genetic files (e.g. VCF, PLINK) is
  COHORT (e.g. `COHORT.vcf.gz`). We will be using this name here, along with 
  several derivatives (e.g. `COHORT_chr1.vcf.gz`).
5. Data have not been subjected to Quality Control (QC).
6. Main QC will be performed with PLINK, some steps may be performed with 
`bcftools`.
7. It is possible that there will be several back and forth format conversions 
during all the process.
8. Data organization (e.g. different directories per step etc.) is the 
responsibility of each partner. We may use a sample structure herein.

## Required software

The tools below will be needed for the proper execution of the analysis plan.
The tools will also be distributed as Docker/Singularity images.

|Tool|Version|Webpage|Download|
|:----|:----|:----|:----|
|bcftools|1.20|[https://www.htslib.org](https://www.htslib.org)|[download](https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2)|
|htslib|1.20|[https://www.htslib.org](https://www.htslib.org)|[download](https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2)|
|VCFtools|0.1.16|[https://vcftools.github.io/](https://vcftools.github.io/)|[download](https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz)|
|PLINK|1.90|[https://www.cog-genomics.org/plink/](https://www.cog-genomics.org/plink/)|[download](https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip)|
|REGENIE|3.5|[https://rgcgithub.github.io/regenie/](https://rgcgithub.github.io/regenie/)|[download](https://github.com/rgcgithub/regenie/releases/download/v3.5/regenie_v3.5.gz_x86_64_Linux.zip)|
|GCTA|1.94.1|[https://yanglab.westlake.edu.cn/software/gcta/](https://yanglab.westlake.edu.cn/software/gcta/)|[download](https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.1-linux-kernel-3-x86_64.zip)|
|MR-MEGA|0.2|[https://genomics.ut.ee/en/tools](https://genomics.ut.ee/en/tools)|[download](https://tools.gi.ut.ee/tools/MR-MEGA_v0.2.zip)|

# 0. Pre-imputation

In the case that genotypic data are not already imputed in HRC version r1.1 or
it is desired that imputation is performed de novo with newer versions of 
imputation/phasing software, the genotypic data should be firstly subjected to
QC. We outline below some steps and assumptions.

## 0.1 Preparation of PLINK files for QC

If you have separate individual ids (IID) and family ids (FID), make sure that 
these are reflected in the names of the samples in VCF and separated by 
underscore (`_`) in the form of `IID_FID`. This means that no underscore should 
be present in the sample name (`IID`) as this will be used for separator by 
PLINK. If you do not have separate family ids, individual ids may contain 
underscores and will be handled appropriately with related `plink` switches. In
**all** cases, rules as described in **General assumption 3** should be 
respected. 

### 0.1.1 Data are in **proper** VCF format

The **proper** VCF specification can be found [here](https://samtools.github.io/hts-specs/VCFv4.2.pdf).
Note that `REF` **must** correspond to the *reference* allele in the genome and 
`ALT` to the *alternative* allele, irrespectively of allele frequencies in the
population. VCF file(s) must also be indexed. If not:

```
tabix COHORT.vcf.gz
```

or, if data are split per chromosome:

```
for CHR in `seq 1 22`
do
  tabix COHORT_chr${CHR}.vcf.gz
done
```

If the data are split per chromosome and you wish to perform the analysis
genome-wide, a genome-wide VCF file should be created in order to proceed with 
PLINK filtering:

```
bcftools concat \
  --output-type z \
  --output COHORT.vcf.gz \
  COHORT_chr{1..22}.vcf.gz
```

Whether you have separate IIDs and FIDs, the  **General assumption 3** should be 
respected. If not, you should create a simple text file with proper names, 
without a header and separated by newline, for example:

Not separate IIDs/FIDs

```
myiid1
myiid2
myiid3
...
myiidn
```

Separate IIDs/FIDs

```
myfid1_myiid1
myfid2_myiid1
myfid3_myiid1
myfid4_myiid2
myfid5_myiid2
...
myfidn_myiidn
```

The file should be called `new_sample_names.txt`. Then, `bcftools reheader` 
should be used to remame the samples:

If you have genome-wide data:

```
bcftools reheader \
  --samples new_sample_names.txt \
  --output COHORT_proper_names.vcf.gz \
  COHORT.vcf.gz
```

If you have data per chromosome:

```
for CHR in `seq 1 22`
do
  bcftools reheader \
    --samples new_sample_names.txt \
    --output COHORT_chr${CHR}_proper_names.vcf.gz \
  COHORT_chr${CHR}.vcf.gz
done
```

If this step has been performed, from now on we assume that `COHORT` is the same
as `COHORT_proper_names` prefix.


If you have proper IIDs and FIDs in the VCF which respect PLINK assumptions
(underscore separation), the command to convert to PLINK format genome-wide is:

```
plink --vcf COHORT.vcf.gz --make-bed --out COHORT
```

or for each chromosome:

```
for CHR in `seq 1 22`
do
  plink --vcf COHORT_chr${CHR}.vcf.gz --make-bed --out COHORT_chr${CHR}
done
``` 

If you don't have separate IIDs and FIDs, then the command for genome-wide data
should be:

```
plink --vcf COHORT.vcf.gz --make-bed --double-id --out COHORT
```

or, per chromosome:

```
for CHR in `seq 1 22`
do
  plink --vcf COHORT_chr${CHR}.vcf.gz --make-bed --double-id \
    --out COHORT_chr${CHR}
done
```

### 0.1.2 Data are in PLINK format

By PLINK format we assume a triplet of BED+BIM+FAM files. If you have the data
in other PLINK supported formats (unlikely) e.g. PED or LGEN, they must be
transformed to BED+BIM+FAM. For example, if you have PED+MAP format:

```
plink --file COHORT --make-bed --out COHORT
```

If you have PLINK files per chromosome and you would prefer (e.g. enough 
computational resources) to perform the analysis genome-wide and not per 
chromosome, these should be merged in one BED+BIM+FAM triplet (FAM should be the 
same for all). It is your responsibility to deal with potential issues that can 
arise from the merging (e.g. duplicates) prior to merging.

```
for CHR in `seq 1 22`
do
  echo COHORT_chr${CHR} >> mergelist.txt
done

plink --merge-list mergelist.txt --make-bed --out COHORT
```

If you prefer to perform the analysis per chromosome (e.g. not enough 
computational resources or enough computational resources to run everything in
parallel), then the above merging step can be skipped. Again, it is your
responsibility to keep track of issues such as duplicates etc.

## 0.2 QC with PLINK

At this point, you should be having a single BED+BIM+FAM triplet

```
COHORT.bed
COHORT.bim
COHORT.fam
```

or a triplet per chromosome (`Z = 1..22`)

```
COHORT_chrZ.bed
COHORT_chrZ.bim
COHORT_chrZ.fam
```

Based on this we proceed with sample and variant filtering.

The following variant filters are recommended:

1. Variant call rate: >98%
2. Minor Allele Frequency: >0.05
3. Hardy-Weinberg equilibrium: >10<sup>-6</sup>

The following sample filters are recommended:

1. Sample call rate: > 95%
2. Heterozygosity: median(heterozygosity) &plusmn; 3 &times; IQR
3. Identity By Descent: > 0.5
4. PCA: outlier removal

Below, we sequentially apply variant and sample filters according to widely
accepted [best practices](https://onlinelibrary.wiley.com/doi/10.1002/sim.6605).

### 0.2.1 Genome-wide QC

First, we calculate heterozygosities:

```
plink \
  --bfile COHORT \
  --out COHORT --het
```

Then, we create a file with sample names with heterozygosities within the limits
of our filter (sample filter #2 above):

```
Rscript \
  -e '{
    hetData <- read.table("COHORT.het",header=TRUE,check.names=FALSE)
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

```
plink \
  --bfile COHORT \
  --out COHORT_tmp \
  --make-bed \
  --geno 0.02 \
  --maf 0.05 \
  --hwe 0.000001 \
  --mind 0.05
```

We now create files with variants and samples to *keep* (samples to keep are 
merged with those passing heterozygosity filters):

```
cut -d" " -f1-2 COHORT_tmp.fam > generic_samples_pass.txt
cut -f2 COHORT_tmp.bim > generic_variants_pass.txt

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
for IBD analysis filtering (can be skipped if not necesseary) and PCA.

```
plink \
  --bfile COHORT \
  --out COHORT_filtered \
  --extract generic_variants_pass.txt \
  --keep all_samples_pass.txt \
  --make-bed
```

In order to perform IBD analysis and filtering (sample filter #3) we firstly 
perform LD-pruning to exclude variants in LD and then IBD calculations with
`plink`. LD-pruning will produce the file `COHORT_filtered.prune.out` which will
be used for IBD calculations:

```
plink \
  --bfile COHORT_filtered \
  --out COHORT_filtered \
  --indep-pairwise 50 5 0.2
  
plink \
  --bfile COHORT_filtered \
  --out COHORT_filtered \
  --genome gz \
  --exclude COHORT_filtered.prune.out
```

The file `COHORT_filtered.genome.gz` is compressed as it may be large. We are
using this in order to find any samples to exclude:

```
Rscript\
  -e '{
    fam <- read.table("COHORT_filtered.fam")
    ibdCoeff <- read.table("COHORT_filtered.genome.gz",header=TRUE)
    ibdCoeff <- ibdCoeff[ibdCoeff$PI_HAT>=0.5,,drop=FALSE]
    if (nrow(ibdCoeff) > 0) {
        bad <- unique(c(ibdCoeff$IID1,ibdCoeff$IID2))
        ii <- match(bad,fam[,2])
        write.table(fam[ii,c(1,2),drop=FALSE],file="ibd_samples_remove.txt",
          col.names=FALSE,row.names=FALSE,quote=FALSE)
    }
  }'
```

Create PLINK files ready for PCA:

```
plink \
  --bfile COHORT_filtered \
  --out COHORT_ibd \
  --remove ibd_samples_remove.txt \
  --make-bed
```

Principal Component Analysis (optional)

*WIP*

### 0.2.2 QC per chromosome and sample

First, we calculate heterozygosities per chromosome:

```
for CHR in `seq 1 22`
do
plink \
  --bfile COHORT_chr${CHR} \
  --out COHORT_chr${CHR} --het
done
```

Then, we create a file with sample names with heterozygosities within the limits
of our filter (sample filter #2 above):

```
Rscript \
  -e '{
    # Read heterozygosity calculations per chromosome
    chrs <- seq_len(22)
    hetData <- lapply(chrs,function(x) {
        hetFile <- paste0("COHORT_chr",x,".het")
        return(read.table(hetFile,header=TRUE,check.names=FALSE))
    })

    # Sum required numbers for each chromosome
    ohom <- rowSums(do.call("cbind",lapply(hetData,function(x) {
        return(x[,3])
    })))
    nnm <- rowSums(do.call("cbind",lapply(hetData,function(x) {
        return(x[,5])
    })))

    # Calculate heterozygosity and filter
    heterozygosity <- 1 - ohom/nnm
    names(heterozygosity) <- rownames(hetData[[1]])
    avg <- median(heterozygosity,na.rm=TRUE)
    dev <- IQR(heterozygosity,na.rm=TRUE)
    keep <- heterozygosity > avg - 3*dev & heterozygosity < avg + 3*dev

    # Write the remaining samples for later
    write.table(hetData[[1]][keep,c("FID","IID"),drop=FALSE],
      file="het_samples_pass.txt",col.names=FALSE,row.names=FALSE,
      quote=FALSE)
  }'
```

The file `het_samples_pass.txt` will be used after the next filters to compile 
a final list of samples to be kept in later analysis.

The following `plink` command will apply variant filters #1,2,3:

```
for CHR in `seq 1 22`
do
plink \
  --bfile COHORT_chr${CHR} \
  --out COHORT_chr${CHR}_tmp \
  --make-bed \
  --geno 0.02 \
  --maf 0.05 \
  --hwe 0.000001
done
```

For sample filter #1, a different approach must be followed, where missingness
is calculated per chromosome and then averaged. If the average values are below 
the inclusion threshold, the sample is excluded. **Please note that the average
values when calculated per chromosome may slightly differ if the analysis is
performed genome-wide**.

Firstly we calculate missing reports for each chromosome:

```
for CHR in `seq 1 22`
do
plink \
  --bfile COHORT_chr${CHR} \
  --out COHORT_chr${CHR} \
  --missing gz
done
```

Then, we will read the report for each chromosome and calculate average 
missingness rates for each sample:

```
Rscript \
  -e '{
    ids <- read.table("COHORT_chr1.imiss.gz",header=TRUE)[,c(1,2)]
    chrs <- seq_len(22)
    missSampleData <- lapply(chrs,function(x) {
        missFile <- paste0("COHORT_chr",x,".imiss.gz")
        missData <- read.table(missFile,header=TRUE)
        fMiss <- missData[,6]
    })
    missingness <- rowMeans(do.call("cbind",missSampleData))
    keep <- missingness < 0.95
    write.table(ids[keep,,drop=FALSE],file="miss_samples_pass.txt",
      col.names=FALSE,row.names=FALSE,quote=FALSE)
  }'
```

The file `miss_samples_pass.txt` contains FID/IIDs of the samples that pass the
missingness test.

We now create files with variants and samples to *keep* (samples to keep are 
merged with those passing heterozygosity filters):

For variants:

```
for CHR in `seq 1 22`
do
  cut -f2 COHORT_chr${CHR}_tmp.bim > generic_variants_pass_chr${CHR}.txt
done
```

For samples:

```
Rscript \
  -e '{
    het <- read.table("het_samples_pass.txt")
    rownames(het) <- het[,2]
    gen <- read.table("miss_samples_pass.txt")
    rownames(gen) <- gen[,2]
    pass <- intersect(rownames(het),rownames(gen))
    write.table(gen[pass,,drop=FALSE],file="all_samples_pass.txt",
      col.names=FALSE,row.names=FALSE,quote=FALSE)
  }'
```

Based on the variant and sample content of the files 
`generic_variants_pass_chrZ.txt` with `Z = 1..22` and `all_samples_pass.txt` 
we create filtered PLINK files for each chromosome. These will be used
for IBD analysis filtering (can be skipped if not necesseary) and PCA.

```
for CHR in `seq 1 22`
do
  plink \
    --bfile COHORT_chr${CHR} \
    --out COHORT_filtered_chr${CHR} \
    --extract generic_variants_pass_chr${CHR}.txt \
    --keep all_samples_pass.txt \
    --make-bed
done
```

In order to perform IBD analysis and filtering (sample filter #3) we firstly 
perform LD-pruning per chromosome to exclude variants in LD and then IBD 
calculations with `plink`. LD-pruning will produce the files 
`COHORT_filtered_chr${CHR}.prune.out` which will be used for IBD calculations:

```
for CHR in `seq 1 22`
do
plink \
  --bfile COHORT_filtered_chr${CHR} \
  --out COHORT_filtered_chr${CHR} \
  --indep-pairwise 50 5 0.2
done

Then we calculate relationships per chromosome:

for CHR in `seq 1 22`
do
plink \
  --bfile COHORT_filtered_chr${CHR} \
  --out COHORT_filtered_chr${CHR} \
  --genome gz \
  --exclude COHORT_filtered_chr${CHR}.prune.out
done
```

The files `COHORT_filtered_chrZ.genome.gz` (`Z = 1..22`) are compressed as they 
may be large. We are using them in order to find any samples to exclude.
**Please note that at present we are not able to provide an official method for
estimating genome-wide IBD coefficients but calculated per chromosome while
using a minimal toolset at the same time for data integrity and interoperability
issues. The following method (median `PH_HAT` accross chromosomes) seems to be
adequately compatible with the original.**

```
Rscript \
  -e '{
    fam <- read.table("COHORT_filtered_chr1.fam")
    chrs <- seq_len(22)
    ibdCoeffs <- do.call("cbind",lapply(chrs,function(x) {
        message("Reading for chromosome ",x)
        ibdFile <- paste0("COHORT_filtered_chr",x,".genome.gz")
        ibdCoeff <- read.table(ibdFile,header=TRUE)
        return(ibdCoeff$PI_HAT)
    }))

    ibdCoeffs <- apply(ibdCoeffs,1,median)
    ibdCoeffNames <- read.table("COHORT_filtered_chr1.genome.gz",header=TRUE)
    ibdCoeffs <- cbind(ibdCoeffNames[,1:4],ibdCoeffs)
    colnames(ibdCoeffs)[5] <- "PI_HAT"

    ibdCoeffs <- ibdCoeffs[ibdCoeffs$PI_HAT>=0.5,,drop=FALSE]
    if (nrow(ibdCoeffs) > 0) {
        bad <- unique(c(ibdCoeffs$IID1,ibdCoeffs$IID2))
        ii <- match(bad,fam[,2])
        write.table(fam[ii,c(1,2),drop=FALSE],file="ibd_samples_remove.txt",
          col.names=FALSE,row.names=FALSE,quote=FALSE)
    }
  }'
```

Create PLINK files ready for PCA:

```
for CHR in `seq 1 22`
do
plink \
  --bfile COHORT_filtered_chr${CHR} \
  --out COHORT_ibd_chr${CHR} \
  --remove ibd_samples_remove.txt \
  --make-bed
done
```

Principal Component Analysis (optional)

*WIP*


## 0.3 Cleanup (optional)

```
rm COHORT_ibd* *prune* COHORT_tmp* 
```

# 1. Post-imputation

The HRC r1.1. 2016 reference panel is well preprocessed, therefore the output
VCF files do not contain variant duplicates and do not contain multi-allelics.
Therefore, their conversion to PLINK after imputation score filtering is more
strightforward. The following assume that the VCF files are named as they are
delivered by the Michigan Imputation Server, i.e. `chr{1..22}.dose.vcf.gz`.

## 1.1 Provision of HRC good quality variants for PCA

You will need to provide to the central analysis team, a list of HRC variant ids
passing imputation quality control. This list will be used along with the lists
from all partners to ensure that the variants used for PCA across all cohorts
will be present across all partners:

```
for CHR in `seq 1 22`
do
  bcftools view \
    --format '%CHROM:%POS' \
    --include 'INFO/R2>0.3' > \
      COHORT_HRC_chr${CHR}_variant.ids &
done
```

The files `COHORT_HRC_chr*_variant.ids` should be provided to the central 
analysis team.

## 1.2 Poorly imputed variant filtering and conversion to PLINK

In the following, we exclude imputed variants with imputation score 
(R<sup>2</sup) less than 0.3. This is done in parallel (per chromosome) where
possible:

```
for CHR in `seq 1 22`
do
  bcftools filter \
    --exclude 'INFO/R2<0.3' \
    --output-type z \
    --output chr${CHR}_filtered.dose.vcf.gz \
    chr${CHR}.dose.vcf.gz &
done
```

Subsequently, we (optionally, depending on the preferred analysis mode, genome-
wide or per chromosome) concatenate the VCF files in preparation for conversion 
to PLINK format for the rest of the filtering as well as PCA operations. We 
assume that the concatenated VCF file is called `COHORT_imputed.vcf.gz`:

```
bcftools concat \
  --output-type z \
  --output COHORT_imputed.vcf.gz \
  chr{1..22}_filtered.dose.vcf.gz
```

Finally, we convert to PLINK:

```
plink --vcf COHORT_imputed.vcf.gz --make-bed --out COHORT_imputed
```

or for each chromosome if we maintain files per chromosome:

```
for CHR in `seq 1 22`
do
  plink --vcf chr${CHR}_filtered.dose.vcf.gz --make-bed \
    --out COHORT_imputed_chr${CHR}
done
```

## 1.3 Post-imputation QC with PLINK

We follow most steps performed in the pre-imputation QC, specifically:

For variant filters:

1. Variant call rate: >98%
2. Minor Allele Frequency: >0.05 (this should also exclude monomorphic variants)
3. Hardy-Weinberg equilibrium: >10<sup>-6</sup>

For sample filters:

1. Sample call rate: > 95%
2. Heterozygosity: median(heterozygosity) &plusmn; 3 &times; IQR
3. Removal of possible gender mismatches
4. PCA: outlier removal

### 1.3.1 Genome-wide QC

As in pre-imputation. we sequentially apply variant and sample filters.

First, we calculate heterozygosities:

```
plink \
  --bfile COHORT_imputed \
  --out COHORT_imputed --het
```

Then, we create a file with sample names with heterozygosities within the limits
of our filter (sample filter #2 above):

```
Rscript \
  -e '{
    hetData <- read.table("COHORT_imputed.het",header=TRUE,check.names=FALSE)
    rownames(hetData) <- hetData$IID
    heterozygosity <- 1 - hetData[,3]/hetData[,5]
    names(heterozygosity) <- rownames(hetData)
    avg <- median(heterozygosity,na.rm=TRUE)
    dev <- IQR(heterozygosity,na.rm=TRUE)
    keep <- heterozygosity > avg - 3*dev & heterozygosity < avg + 3*dev
    goodHet <- rownames(hetData)[keep]
    write.table(hetData[keep,c("FID","IID"),drop=FALSE],
      file="het_samples_ipass.txt",col.names=FALSE,row.names=FALSE,quote=FALSE)
  }'
```

The file `het_samples_ipass.txt` will be used after the next filters to compile 
a final list of samples to be kept in later analysis.

The following `plink` command will apply variant filters #1,2,3 and sample 
filter #1:

```
plink \
  --bfile COHORT_imputed \
  --out COHORT_imputed_tmp \
  --make-bed \
  --geno 0.02 \
  --maf 0.05 \
  --hwe 0.000001 \
  --mind 0.05
```

We now create files with variants and samples to *keep* (samples to keep are 
merged with those passing heterozygosity filters):

```
cut -d" " -f1-2 COHORT_imputed_tmp.fam > generic_samples_ipass.txt
cut -f2 COHORT_imputed_tmp.bim > generic_variants_ipass.txt

Rscript \
  -e '{
    het <- read.table("het_samples_ipass.txt")
    rownames(het) <- het[,2]
    gen <- read.table("generic_samples_ipass.txt")
    rownames(gen) <- gen[,2]
    pass <- intersect(rownames(het),rownames(gen))
    write.table(gen[pass,,drop=FALSE],file="all_samples_ipass.txt",
      col.names=FALSE,row.names=FALSE,quote=FALSE)
  }'
```

Create PLINK files ready for PCA:

```
plink \
  --bfile COHORT_imputed \
  --out COHORT_imputed_filtered \
  --extract generic_variants_ipass.txt \
  --keep all_samples_ipass.txt \
  --make-bed
```


Principal Component Analysis

*WIP*

### 1.3.2 QC per chromosome and sample

First, we calculate heterozygosities per chromosome:

```
for CHR in `seq 1 22`
do
plink \
  --bfile COHORT_imputed_chr${CHR} \
  --out COHORT_imputed_chr${CHR} --het
done
```

Then, we create a file with sample names with heterozygosities within the limits
of our filter (sample filter #2 above):

```
Rscript \
  -e '{
    # Read heterozygosity calculations per chromosome
    chrs <- seq_len(22)
    hetData <- lapply(chrs,function(x) {
        hetFile <- paste0("COHORT_imputed_chr",x,".het")
        return(read.table(hetFile,header=TRUE,check.names=FALSE))
    })

    # Sum required numbers for each chromosome
    ohom <- rowSums(do.call("cbind",lapply(hetData,function(x) {
        return(x[,3])
    })))
    nnm <- rowSums(do.call("cbind",lapply(hetData,function(x) {
        return(x[,5])
    })))

    # Calculate heterozygosity and filter
    heterozygosity <- 1 - ohom/nnm
    names(heterozygosity) <- rownames(hetData[[1]])
    avg <- median(heterozygosity,na.rm=TRUE)
    dev <- IQR(heterozygosity,na.rm=TRUE)
    keep <- heterozygosity > avg - 3*dev & heterozygosity < avg + 3*dev

    # Write the remaining samples for later
    write.table(hetData[[1]][keep,c("FID","IID"),drop=FALSE],
      file="het_samples_pass.txt",col.names=FALSE,row.names=FALSE,
      quote=FALSE)
  }'
```

The file `het_samples_pass.txt` will be used after the next filters to compile 
a final list of samples to be kept in later analysis.

The following `plink` command will apply variant filters #1,2,3:

```
for CHR in `seq 1 22`
do
plink \
  --bfile COHORT_imputed_chr${CHR} \
  --out COHORT_imputed_chr${CHR}_tmp \
  --make-bed \
  --geno 0.02 \
  --maf 0.05 \
  --hwe 0.000001
done
```

For sample filter #1, a different approach must be followed, where missingness
is calculated per chromosome and then averaged. If the average values are below 
the inclusion threshold, the sample is excluded. **Please note that the average
values when calculated per chromosome may slightly differ if the analysis is
performed genome-wide**.

Firstly we calculate missing reports for each chromosome:

```
for CHR in `seq 1 22`
do
plink \
  --bfile COHORT_imputed_chr${CHR} \
  --out COHORT_imputed_chr${CHR} \
  --missing gz
done
```

Then, we will read the report for each chromosome and calculate average 
missingness rates for each sample:

```
Rscript \
  -e '{
    ids <- read.table("COHORT_imputed_chr1.imiss.gz",header=TRUE)[,c(1,2)]
    chrs <- seq_len(22)
    missSampleData <- lapply(chrs,function(x) {
        missFile <- paste0("COHORT_imputed_chr",x,".imiss.gz")
        missData <- read.table(missFile,header=TRUE)
        fMiss <- missData[,6]
    })
    missingness <- rowMeans(do.call("cbind",missSampleData))
    keep <- missingness < 0.95
    write.table(ids[keep,,drop=FALSE],file="miss_samples_pass.txt",
      col.names=FALSE,row.names=FALSE,quote=FALSE)
  }'
```

The file `miss_samples_pass.txt` contains FID/IIDs of the samples that pass the
missingness test.

We now create files with variants and samples to *keep* (samples to keep are 
merged with those passing heterozygosity filters):

For variants:

```
for CHR in `seq 1 22`
do
  cut -f2 COHORT_imputed_chr${CHR}_tmp.bim > generic_variants_pass_chr${CHR}.txt
done
```

For samples:

```
Rscript \
  -e '{
    het <- read.table("het_samples_pass.txt")
    rownames(het) <- het[,2]
    gen <- read.table("miss_samples_pass.txt")
    rownames(gen) <- gen[,2]
    pass <- intersect(rownames(het),rownames(gen))
    write.table(gen[pass,,drop=FALSE],file="all_samples_pass.txt",
      col.names=FALSE,row.names=FALSE,quote=FALSE)
  }'
```

Based on the variant and sample content of the files 
`generic_variants_pass_chrZ.txt` with `Z = 1..22` and `all_samples_pass.txt` 
we create filtered PLINK files for each chromosome. These will be used
for IBD analysis filtering (can be skipped if not necesseary) and PCA.

```
for CHR in `seq 1 22`
do
  plink \
    --bfile COHORT_imputed_chr${CHR} \
    --out COHORT_imputed_filtered_chr${CHR} \
    --extract generic_variants_pass_chr${CHR}.txt \
    --keep all_samples_pass.txt \
    --make-bed
done
```

In order to perform IBD analysis and filtering (sample filter #3) we firstly 
perform LD-pruning per chromosome to exclude variants in LD and then IBD 
calculations with `plink`. LD-pruning will produce the files 
`COHORT_imputed_filtered_chr${CHR}.prune.out` which will be used for IBD
calculations:

```
for CHR in `seq 1 22`
do
plink \
  --bfile COHORT_imputed_filtered_chr${CHR} \
  --out COHORT_imputed_filtered_chr${CHR} \
  --indep-pairwise 50 5 0.2
done

Then we calculate relationships per chromosome:

for CHR in `seq 1 22`
do
plink \
  --bfile COHORT_imputed_filtered_chr${CHR} \
  --out COHORT_imputed_filtered_chr${CHR} \
  --genome gz \
  --exclude COHORT_imputed_filtered_chr${CHR}.prune.out
done
```

The files `COHORT_filtered_chrZ.genome.gz` (`Z = 1..22`) are compressed as they 
may be large. We are using them in order to find any samples to exclude.
**Please note that at present we are not able to provide an official method for
estimating genome-wide IBD coefficients but calculated per chromosome while
using a minimal toolset at the same time for data integrity and interoperability
issues. The following method (median `PH_HAT` accross chromosomes) seems to be
adequately compatible with the original.** Also, please note that for large
populations, reading IBD `*.genome.gz` files may take a considerable amount of
time. Consider replacing `read.table` with `fread` from the package 
`data.table`.

```
Rscript \
  -e '{
    fam <- read.table("COHORT_imputed_filtered_chr1.fam")
    chrs <- seq_len(22)
    ibdCoeffs <- do.call("cbind",lapply(chrs,function(x) {
        message("Reading for chromosome ",x)
        ibdFile <- paste0("COHORT_imputed_filtered_chr",x,".genome.gz")
        ibdCoeff <- read.table(ibdFile,header=TRUE)
        return(ibdCoeff$PI_HAT)
    }))

    ibdCoeffs <- apply(ibdCoeffs,1,median)
    ibdCoeffNames <- read.table("COHORT_imputed_filtered_chr1.genome.gz",
        header=TRUE)
    ibdCoeffs <- cbind(ibdCoeffNames[,1:4],ibdCoeffs)
    colnames(ibdCoeffs)[5] <- "PI_HAT"

    ibdCoeffs <- ibdCoeffs[ibdCoeffs$PI_HAT>=0.5,,drop=FALSE]
    if (nrow(ibdCoeffs) > 0) {
        bad <- unique(c(ibdCoeffs$IID1,ibdCoeffs$IID2))
        ii <- match(bad,fam[,2])
        write.table(fam[ii,c(1,2),drop=FALSE],file="ibd_samples_remove.txt",
          col.names=FALSE,row.names=FALSE,quote=FALSE)
    }
  }'
```

Create reduced (by removed samples) PLINK files ready for PCA:

```
for CHR in `seq 1 22`
do
plink \
  --bfile COHORT_imputed_filtered_chr${CHR} \
  --out COHORT_imputed_ibd_chr${CHR} \
  --remove ibd_samples_remove.txt \
  --make-bed
done
```

And then merge them in preparation for PCA and/or PC projection to 1000 genomes:

```
for CHR in `seq 1 22`
do
  echo COHORT_imputed_ibd_chr${CHR} >> mergelist.txt
done

plink --merge-list mergelist.txt --out COHORT_imputed_filtered_merged
```

Principal Component Analysis

*WIP*

## 1.4 Non-genetic related sample exclusions

Individuals of all age groups will be included in the analyses. However, the
following samples should also be excluded from any analyses:

- Outliers (based on 5 &times; SD) for any trait will be excluded.
- Known Mendelian cases of obesity will be excluded.
- Pregnant women and individuals with any other acute or chronic clinical 
condition that significantly alters normal body weight will be excluded.

**IMPORTANT**: Analyses will take place separately for children/adolescents <18 
years old and adults older than 18 years old.

## 1.5 Principal Component projections

The central analysis team will provide three files based on Principal Component
Analysis of the current 1000 genomes data. These files will be:

- `pca_variants.txt`: a set of SNPs whose loadings will be used to project your
local genotypes to 1000 genomes and be used as covariates in subsequent GWAS.
- `loads_1000g.txt`: the aforementioned loadings.
- `means_1000g.txt`: means and standard deviations of the loadings required by
`flashpca`.

In this step the file `pca_variants.txt` will be used in order to create a PLINK
fileset with a subset of variants contained in this file. Then, these will be
used with the files `loads_1000g.txt` and `means_1000g.txt` along with the 
projection functionality of `flashpca` to create the PC covariates for your
cohort. The PC projections have to be provided to the central analysis team.

**IMPORTANT**: Although the central analysis team will make sure that the 
variants  used for PCA projections are common across all partners (see step 1.1
above), it remains **your** responsibility to make sure that **all** the 
variants in the provided file `pca_variants.txt` are present in your cohort. If
the cohort is not very small and has been imputed with HRC panel, this should
not be a problem. Otherwise, please contact the central analysis team 
**immediately**.

Firstly, create the PLINK files for PCA projection:

```
plink \
  --bfile COHORT_imputed_filtered_merged \
  --out COHORT_for_PCA \
  --extract pca_variants.txt \
  --make-bed
```

Then, project:

```
flashpca \
  --bfile COHORT_for_PCA \
  --inmeansd means_1000g.txt \
  --inload loads_1000g.txt
  --project \
  --outproj projections.txt \
  --verbose
```

We are now editing the file `projections.txt` so as to:

1. Anonymize the individual ids accompanying the PC projections
2. Create a map between the new and the old ids

```
Rscript \
  -e '{
    NUMBERS <- as.character(seq(0,9))

    .randomString <- function(n=1,s=5) {
        SPACE <- sample(c(LETTERS,NUMBERS),length(LETTERS)+length(NUMBERS))
        return(do.call(paste0,replicate(s,sample(SPACE,n,replace=TRUE),FALSE)))
    }

    PREFIX <- "YOUR_PARTNER_ID_IN_GRANT e.g. HUA"
    PREFIX <- "HUA"

    # Read in projections
    proj <- read.delim("projections.txt")

    # Strip FIDs
    proj <- proj[,-1,drop=FALSE]

    # Create random names and the map
    # We create 10-letter random strings to allow enough variability and avoid
    # random same string in larger populations
    set.seed(42)
    newids <- paste(PREFIX,.randomString(nrow(proj),10),sep="")
    map <- cbind(proj$IID,newids)
    colnames(map) <- c("cohort_iid","pseudo_iid")

    # Replace the ids and write output
    proj$IID <- newids
    write.table(map,file="cohort_pseudo_map.txt",row.names=FALSE,sep="\t",
        quote=FALSE)
    write.table(proj,file="shared_projections.txt",row.names=FALSE,sep="\t",
        quote=FALSE)

  }'
```

The file `shared_projections.txt` will be returned to the central analysis team 
for processing. Its totality or part of it (e.g. if samples will be exluded as 
outliers) will be returned to you to be used along with other covariates in the
subsequent GWAS.

The file `cohort_pseudo_map.txt` must be safely stored by you to map the 
included individuals returned by the central analysis team (if anyone excluded)
for the subsequent site-specific GWAS. The `set.seed()` function ensures that
the same ids can be regenerated in case of loss.

# 2. Genome-wide association analyses

## 2.1 General

Please use the XXX script to perform the genome-wide association analyses 
required herein. 

* An additive model of inheritance will be assumed.
* In case you have a case-control study, and the disease phenotype is associated 
with anthropometric data, please make sure you run the analysis separately for 
cases and controls.
* In case of a family study or cryptic relatedness, please include a kinship 
matrix in your analysis using XXX.
* Any other cohort-specific covariates should be included as covariates in the 
models.

The following analyses will take place:

1.

$$ BMI(t) = age(t) + age(t)^2 + sex + \sum_{i=1}^n PC_i + OtherCohortSpecificCovariates $$

2.

$$ \Delta BMI(t) = \Delta age(t) + sex + \sum_{i=1}^n PC_i + OtherCohortSpecificCovariates $$

3.

$$ BodyWeight(t) = age(t) + age(t)^2 + sex + \sum_{i=1}^n PC_i + OtherCohortSpecificCovariates $$

4. Others *WIP*

## 2.2 Analysis with REGENIE

## 2.2.1 Preparation of covariate and phenotype files

The covariate files should include covariates for each analysis mentioned above
plus the PCs. Based on the accompanying toy dataset from HUA, we provide an
example for case (i).

For the covariate file, we need the phenotypes as well as the PCs. The PCs will
be provided by the central analysis team. For this example, we calculate 10 PCs
with PLINK.

Firstly, perform pruning (will be done based on the reference panel for the
canonical analysis):

```
plink \
  --bfile toy \
  --out toy \
  --indep-pairwise 50 5 0.2
```

Then, calculate 10 PCs exluding the pruned SNPs:

```
plink \
  --bfile toy \
  --out toy \
  --pca 10 header \
  --exclude toy.prune.out
```

The output is written to `toy.eigenvec`. This is just for example purposes. 
During the normal execution, you will not calculate PCs with PLINK but rather 
use the PC projections file provided by the central analysis team. We then
construct the covariates and phenotype file:

```
Rscript \
  -e '{
    # Read all phenotypes
    phen <- read.delim("toy.txt")
    
    # Construct 1st part of covariates files
    rownames(phen) <- phen$iid
    phen$age2 <- phen$age^2
    covs <- data.frame(
      FID=phen$fid,
      IID=phen$iid,
      sex=phen$sex,
      age=phen$age,
      age2=phen$age2
    )
    
    # 2nd part - PCs
    pcs <- read.table("toy.eigenvec",header=TRUE)
    rownames(pcs) <- pcs$IID
    pcs <- pcs[rownames(phen),]
    pcs <- pcs[,grep("PC",colnames(pcs))]
    covs <- cbind(covs,pcs)
    
    # Write the covariates file
    write.table(covs,file="toy_regenie_covariates.txt",row.names=FALSE,
        quote=FALSE)
    
    # Construct the phenotype file
    resp <- data.frame(
      FID=phen$fid,
      IID=phen$iid,
      bmi=phen$bmi_b
    )
    write.table(resp,file="toy_regenie_phenotype.txt",row.names=FALSE,
        quote=FALSE)
  }'
```

**NOTE**: We just put the `bmi_b` variable in the phenotype file, however, we
could have put all the variables of interest for multi-trait analysis since
`regenie` supports it.

## 2.2.2 Execute REGENIE with the toy dataset

Step 0: exclude monomorphic and low variance SNPs from the toy dataset as they
cause `regenie` to crash. Generally, QC **must** be performed prior to running
`regenie`:

```
plink \
  --bfile toy \
  --out toyf \
  --make-bed \
  --maf 0.01
```

Step 1:

```
regenie \
  --step 1 \
  --bed toyf \
  --covarFile toy_regenie_covariates.txt \
  --phenoFile toy_regenie_phenotype.txt \
  --bsize 100 \
  --out fit_bmi_out
```

Step 2:

```
regenie \
  --step 2 \
  --bed toyf \
  --covarFile toy_regenie_covariates.txt \
  --phenoFile toy_regenie_phenotype.txt \
  --bsize 200 \
  --pThresh 0.05 \
  --pred fit_bmi_out_pred.list \
  --ignore-pred \
  --out test_bmi_out_firth
```

The summary statistics are in `test_bmi_out_firth_bmi.regenie`.

**NOTE**: 

1. All the parameters will be present in the final analysis plan and decided by
the central analysis team.
2. We used `--ignore-pred` here because our sample for Step 1 is too small to 
produce whole genome predictions required by `regenie`. It should **not** be the
case with real data.

## 2.3 Analysis with GCTA

*WIP*

## Notes

* If you work in a multicore Linux system with enough available RAM, use the
ambersand (`&`) symbol at the end of command lines, for example:

```
for CHR in `seq 1 22`
do
plink \
  --bfile COHORT_chr${CHR} \
  --out COHORT_chr${CHR} --het &
done
```

This will spawn 22 processes working in parallel and will produce one 
heterozygosity file per chromosome.

* Consider using `nohup` for longer calculations (e.g. many samples on imputed
data) that cannot finish in an interactive shell session in reasonable time. For
example, put the commands to run in a file `commands.sh`, for example:

```
#!/bin/bash

for CHR in `seq 1 22`
do
plink \
  --bfile COHORT_chr${CHR} \
  --out COHORT_chr${CHR} --het &
done
```

Then, execute:

```
nohup sh commands.sh > commands.log &
```

and you can safely close the session and check e.g. the next day.

## Open issues

* Which program will be used for GWAS? We have drafted examples using REGENIE.
It was suggested by R. Pool to use GCTA as SAIGE and REGENIE do not behave well
with twins.
* "Anonymization" of PCs should be checked/verified.
* Do we need also `FID` for the clustering of projected PCs? If yes the process
above must be amended.
* Do we need other population information to be provided along with the PC
projections?
* Will [the inverse normal transformation of residuals](https://www.biostars.org/p/312945/) 
will be used for any BMI related associations?
* Other intermediate steps?
