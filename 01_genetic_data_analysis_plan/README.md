Draft BETTER4U analysis plan
================================================================================

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

* [bcftools 1.20](https://samtools.github.io/bcftools/) - [download](https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2) - [manual](https://samtools.github.io/bcftools/bcftools.html)
* [htslib 1.20](https://www.htslib.org/) - [download](https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2) - [manual](https://www.htslib.org/doc/bgzip.html)
* [PLINK 1.90](https://www.cog-genomics.org/plink/) - [download](https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip) - [manual](https://www.cog-genomics.org/plink/1.9/index)

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
Note that `REF` **must** correspond to the *reference* allele in the genoe and 
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

If the data are split per chromosome, a genome-wide VCF file should be created
in order to proceed with PLINK filtering:

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
myiid1_myfid1
myiid2_myfid1
myiid3_myfid1
myiid4_myfid2
myiid5_myfid2
...
myiidn_myfidn
```

The file should be called `new_sample_names.txt`. Then, `bcftools reheader` 
should be used to remame the samples:

```
bcftools reheader \
  --samples new_sample_names.txt \
  --output COHORT_proper_names.vcf.gz \
  COHORT.vcf.gz
```

If you have proper IIDs and FIDs in the VCF which respect PLINK assumptions
(underscore separation), the command to convert to PLINK format is:

```
plink --vcf COHORT.vcf.gz --make-bed --out COHORT
```

If you don't have separate IIDs and FIDs, then the command should be:

```
plink --vcf COHORT.vcf.gz --make-bed --double-id --out COHORT
```

### 0.1.2 Data are in PLINK format

By PLINK format we assume a triplet of BED+BIM+FAM files. If you have the data
in other PLINK supported formats (unlikely) e.g. PED or LGEN, they must be
transformed to BED+BIM+FAM. For example, if you have PED+MAP format:

```
plink --file COHORT --make-bed --out COHORT
```

If you have PLINK files per chromosome, these should be merged in one 
BED+BIM+FAM triplet (FAM should be the same for all). It is your responsibility
to deal with potential issues that can arise from the merging (e.g. duplicates)
prior to merging.

```
for CHR in 1..22
do
  echo COHORT_chr${CHR} >> mergelist.txt
done

plink --merge-list mergelist.txt --make-bed --out COHORT
```

## 0.2 QC with PLINK

At this point, you should be having a single BED+BIM+FAM triplet:

```
COHORT.bed
COHORT.bim
COHORT.fam
```

Based on this we proceed with sample and variant filtering.

The following variant filters are recommended:

1. Variant call rate: >98%
2. Minor Allele Frequency: >0.05
3. Hardy-Weinberg equilibrium: >10<sup>-6</sup>

The following sample filters are recommended:

1. Sample call rate: > 95%
2. Heterozygosity: median(heterozygosity) +/- 3IQR
3. Identity By Descent: > 0.5
4. PCA: outlier removal

Below, we sequentially apply variant and sample filters according to widely
accepted [best practices](https://onlinelibrary.wiley.com/doi/10.1002/sim.6605).

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

We now create files with variants and samples to *keep* (samples to keep will 
are merged with those passing heterozygosity filters):

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

Principal Component Analysis

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

## 1.1 Poorly imputed variant filtering and conversion to PLINK

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

Subsequently, we concatenate the VCF files in preparation for conversion to
PLINK format for the rest of the filtering as well as PCA operations. We assume
that the concatenated VCF file is called `COHORT_imputed.vcf.gz`:

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

## 1.2 Post-imputation QC with PLINK

We follow most steps performed in the pre-imputation QC, specifically:

For variant filters:

1. Variant call rate: >98%
2. Minor Allele Frequency: >0.05 (this should also exclude monomorphic variants)
3. Hardy-Weinberg equilibrium: >10<sup>-6</sup>

For sample filters:

1. Sample call rate: > 95%
2. Heterozygosity: median(heterozygosity) +/- 3 x IQR
3. Removal of possible gender mismatches
4. PCA: outlier removal

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


## Notes

If multicore, use `&`
Consider using `nohup`

