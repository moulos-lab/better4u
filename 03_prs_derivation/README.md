BETTER4U Task 3.7 analysis plan (internal and shared)
================================================================================

## Authors

**Panagiotis Moulos** <sup>1,2</sup><br/>
<sup>1</sup>Harokopio University of Athens, <sup>2</sup>BSRC Alexander Fleming

### Contact details

* Panagiotis Moulos (moulos@fleming.gr)
* Jon Anders Eriksson (anders.eriksson@ut.ee)
* Nana Kalafati (nkalafati@gmail.com)

# Introduction

This analysis plan covers the work required for the completion of BETTER4U WP3
Task 3.7: Development of an iterative PRS extraction pipeline

The task comprises X main parts (from PPT)

## Required software

* [PLINK 1.90](https://www.cog-genomics.org/plink/)
* [samtools](https://github.com/samtools/samtools/releases/download/1.22/samtools-1.22.tar.bz2)
* [bcftools](https://samtools.github.io/bcftools/bcftools.html)
* [SnpEff/SnpSift](https://pcingola.github.io/SnpEff/download/)
* [GCTB](https://cnsgenomics.com/software/gctb/#Download)

## Installation of SBayesRC, PRS-CS and resources

### Directory structure

The proposed work assumes the following directory structure:

```
WORKSPACE
|__resources
|  |__LD
|  |_...
|__work
|  |__METAL
|__|__PRS
|__bootstrap
|  |__1
|__|__2
|__|__...
```

We create the structure:

```bash
export WORKSPACE=/media/storage3/playground/b4uprs

mkdir -p $WORKSPACE/resources
mkdir -p $WORKSPACE/work/METAL
mkdir -p $WORKSPACE/work/PRS/baseline
mkdir -p $WORKSPACE/work/PRS/bootstrap
mkdir -p $WORKSPACE/work/LOO
```

### SBayesRC

We follow the instructions at the [SBayesRC](Follow instructions at https://github.com/zhilizheng/SBayesRC) page. A known bug requires downgrading
(temporarily) the BOOST library headers for R. Also, parallelization does not
work well in the vanilla version. We created a [fork](https://github.com/pmoulos/SBayesRC)
that uses the R package `parallel` inspired by another [fork](https://github.com/andrew-terpolovsky/SBayesRC).

```r
library(remotes)

# Downgrade boost headers...
install.packages("https://cran.r-project.org/src/contrib/Archive/BH/BH_1.81.0-1.tar.gz")
install.packages(c("Rcpp","data.table","stringi","BH","RcppEigen"))
#remotes::install_github("zhilizheng/SBayesRC")
remotes::install_github("pmoulos/SBayesRC")
```

### PRS-CS

We follow the instructions at the [PRS-CS](https://github.com/getian107/PRScs)
page.

```bash
pip install h5py scipy

cd $WORKSPACE/resources

git clone https://github.com/getian107/PRScsx.git
```

### Resources

#### PRS-CS resources

We follow the instructions at the [PRS-CS](https://github.com/getian107/PRScs)
page.

```bash
cd $WORKSPACE/resources

mkdir PRScsLD && cd PRScsLD

wget https://www.dropbox.com/s/mt6var0z96vb6fv/ldblk_1kg_eur.tar.gz?dl=0 -O ldblk_1kg_eur.tar.gz
wget https://www.dropbox.com/s/rhi806sstvppzzz/snpinfo_mult_1kg_hm3?dl=0 -O snpinfo_mult_1kg_hm3
tar -xvf ldblk_1kg_eur.tar.gz
cd ..
```

#### SBayesRC resources

The following resources are required by SBayesRC. The sbayes server is slow, 
therefore we suggest downloading the resources from their Google Drive using
`gdown`.

```bash
cd $WORKSPACE/resources

gdown 1mKQ3uU_XD6zlNefxEWMl1M42I0gTyOs3 -O ukbEUR_Imputed.tar.xz
gdown 1-dUPvduYB1zZewsItCNKcM7RvOAOojlP -O annot_baseline2.2.zip

tar -xvf ukbEUR_Imputed.zip # Takes long!
unzip annot_baseline2.2.zip
```

#### 1000 genomes VCF files

[VCF files per chromosomes](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502) 
from the latest 1000 genomes project FTP (phase 3)

# Preparation of LD structure from 1000 genomes for usage with SBayesRC

## Process

### 1. Download the files as VCF.gz (and tab-indices) and PED file

```bash
cd $WORKSPACE/resources
mkdir TGP && cd TGP

PREFIX="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr"
SUFFIX=".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

for CHR in `seq 1 22`; 
do
    wget ${PREFIX}${CHR}${SUFFIX} ${PREFIX}${CHR}${SUFFIX}".tbi";
done

# PED file to extract EUR samples
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped
```

### 2. Download the GRCh37/hg19 reference genome

For better data integrity, we download the human reference genome version used
by the consortium:

```
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
gunzip human_g1k_v37.fasta.gz
samtools faidx human_g1k_v37.fasta
```

### 3. Extract EUR samples from PED file to be used for filtering VCFs

```bash
# Get the samples from PED
grep -P 'CEU|FIN|GBR|IBS|TSI' 20130606_g1k.ped | \
  perl -lane '$seen{$F[1]}++ || print $F[1]' > ped_eur.ids
  
# Some do not exist in the VCF files...
bcftools query -l \
  ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz > \
  all_samples.ids

Rscript \
  -e '{
    a <- read.table("ped_eur.ids")[,1]
    b <- read.table("all_samples.ids")[,1]
    ab <- intersect(a,b)
    writeLines(ab,"eur.ids")
  }'
```

### 4. Convert the 1000 Genomes files to BCF

In the following we convert the VCF files from 1000 genomes data to BCF files
which is suitable for further dowstream analysis and filtering prior to LD
construction. Specifically:

* Include only bi-allelic SNPs with MAF > 0.01 (1st pipe)
* Inlcude only samples of European ancestry (1st pipe)
* Ensure that multi-allelic calls are split and that indels are left-aligned 
compared to reference genome (2nd pipe)
* Remove any variant IDs to later be re-annotated (3rd pipe)
* Remove duplicates (4th pipe)

```
for CHR in `seq 1 22`; 
do
  echo "Converting $CHR"
  
  bcftools view --include 'INFO/AF > 0.01' --types snps \
    --min-alleles 2 --max-alleles 2 --samples-file eur.ids \
    ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz | \
  bcftools norm --multiallelics -any --check-ref w \
    --fasta-ref human_g1k_v37.fasta | \
  bcftools annotate --remove ID | \
  bcftools norm --output-type z --rm-dup both \
    --output ALL_1000G_chr${CHR}.vcf.gz &
done

wait

for CHR in `seq 1 22`; 
do
  tabix ALL_1000G_chr${CHR}.vcf.gz &
done
```

### 5. Annotate with SnpEff and dbSNP for hg19

Annotate with SnpSift in parallel...

```bash
for CHR in `seq 1 22`; 
do
  java -jar $WORKSPACE/resources/snpEff/SnpSift.jar annotate -noInfo -noLog \
    -id $WORKSPACE/resources/dbsnp/dbSNP155.vcf.gz \
    ALL_1000G_chr${CHR}.vcf.gz |
  bcftools view --include 'ID != "."' --output-type z \
    --output ALL_1000G_chr${CHR}_dbsnp_nomiss.vcf.gz &
done
```

...and compress the output VCFs

```bash
for CHR in `seq 1 22`; 
do
  bgzip ALL_1000G_chr${CHR}_dbsnp_nomiss.vcf &
done
```

### 6. Merge all VCFs to one file for later LD

```bash
NEWPRE="ALL_1000G_chr"
NEWSUF="_dbsnp_nomiss.vcf.gz"
bcftools concat --output-type z \
  --output ALL_1000G_dbsnp_nomiss.vcf.gz ${NEWPRE}{1..22}${NEWSUF}
tabix ALL_1000G_dbsnp_nomiss.vcf.gz
```

### 7. Convert the annotated file to PLINK for downstream LD panel

```bash
plink \
  --vcf ALL_1000G_dbsnp_nomiss.vcf.gz \
  --keep-allele-order \
  --vcf-idspace-to _ \
  --const-fid \
  --allow-extra-chr 0 \
  --split-x b37 no-fail \
  --make-bed \
  --out ALL_1000G_dbsnp_nomiss
```

### 8. Use SBayesRC to generate EUR 1000 genomes LD

SBayesRC script uses GCTB software to generate LD blocks, so GCTB should be
installed and accessible from the main path:

```r
Rscript \
  -e '{
    library(SBayesRC)
    library(parallel)

    WORKSPACE <- Sys.getenv("WORKSPACE")

    maFile <- file.path(WORKSPACE,"work/METAL/metal_b4u.ma")
    genotype <- file.path(WORKSPACE,"resources/TGP/ALL_1000G_dbsnp_nomiss")
    outDir <- file.path(WORKSPACE,"resources/EUR_LD")
    threads <- 16

    # Step 1
    SBayesRC::LDstep1(
        mafile=maFile,
        genoPrefix=genotype,
        outDir=outDir
    )
  }'
```

A current [bug](https://github.com/zhilizheng/SBayesRC/issues/55) in SBayesRC LD
constructing wrapper is causing simultaneous spawn of 591(!) processes leading
to memory crash, so we need to fix.

```bash
sed -i 's/&> \(\S*\.log\)/> \1 2>\&1/' $WORKSPACE/resources/EUR_LD/ld.sh
```

Back to R now:

```r
Rscript \
  -e '{
    library(SBayesRC)
    library(parallel)

    WORKSPACE <- Sys.getenv("WORKSPACE")

    maFile <- file.path(WORKSPACE,"work/METAL/metal_b4u.ma")
    genotype <- file.path(WORKSPACE,"resources/TGP/ALL_1000G_dbsnp_nomiss")
    outDir <- file.path(WORKSPACE,"resources/EUR_LD")
    threads <- 16

    # Step 2
    mclapply(seq(591),function(i) {
        SBayesRC::LDstep2(
            outDir=outDir,
            blockIndex=i
        )
    },mc.cores=threads)

    # Step 3
    mclapply(seq(591),function(i) {
        SBayesRC::LDstep3(
            outDir=outDir,
            blockIndex=i
        )
    },mc.cores=threads)

    # Step 4
    SBayesRC::LDstep4(outDir=outDir)
  }'
```

### 9. Cleaning

```bash
rm -r $WORKSPACE/resources/EUR_LD/*.log $WORKSPACE/resources/EUR_LD/*.full.bin \
  $WORKSPACE/resources/EUR_LD/*.full.info $WORKSPACE/resources/EUR_LD/snplist \
  $WORKSPACE/resources/EUR_LD/ld.sh
```

# PRS processes

The following describe the derivation and evaluation process for baseline PRS 
for BMI and weight/BMI change derivation. Bootstrapping stabilization will be
performed only for the best baseline PRS.

## Baseline PRS derivation with SBayesRC and PRS-CS

We assume that we have one file of METAL derived summary statistics for each 
chromosome called `metal.chrZ.txt` with Z between 1 and 22. These are stored in
`$WORKSPACE/work/METAL`. 

### Data preparation

Based on the [README](https://docs.google.com/document/d/1wolKXW5SfeNL4RtxCfcGj6JmC992r2cJ8b_Vw4ge2OE/edit?usp=drive_link) 
of the METAL analysis performed by Anders Eriksson, there were several data
transformations which require further pre-processing prior to running any of the
tools for PRS derivation. It is stated that:

>GWAS summary statistics files were harmonised with respect to reference and alternate allele, so that effect allele is the alternate allele.
>SNPs were subset to MTAG SNPs for downstream multi-trait analyses.
>METAL analyses was carried out to combine summary statistics from the different BETTER4U partners.
>First, P values were adjusted for inflation (fitting a second degree polynomial of expected vs. observed P values for P > 0.01 and using the linear component to scale P values). Second, meta-analyses were carried out using (inverse variance weighted) METAL: new Z values were calculated using inverse transform of the adjusted P-values and the sign of the effect size. Z scores were combined using the sample size as weights.

Therefore:

* p-values were adjusted for inflation before meta-analysis
* Then, Z were recomputed from those adjusted p-values, keeping the sign of the 
effect.
* METAL was run in a sample-size weighted mode, not truly inverse-variance
* The `SE` column is a pseudo-SE and can be negative, because it’s not really a 
standard error at all.

As a result of the points above, betas and SEs must be reconstructed from Z, N
and AF. The reconstructed betas and SEs will reflect **standardized phenotype**
scales (inerpretation per SD unit of the phenotype per additional ALT allele) 
which is inline with how GCTA `--fastGWA` works internally and also how PRS 
tools expect effects (per allele and not per standardized genotype). 

Specfically, the following conversions are made, assuming an additive model and 
small per-SNP R<sup>2</sup>, with genotypes coded 0/1/2 and effect allele=ALT,
let:

$$
\text{Var}(G)=2f(1-f), \quad f=\text{ALT allele frequency (AF)}, \quad N=\text{per-SNP (effective) sample size}
$$

then

$$
SE \approx\ \frac{1}{\sqrt{N\text{Var}(G)}} \quad\text{and}\quad
\hat\beta \approx \frac{Z}{\sqrt{N\text{Var}(G)}} = Z \cdot SE
$$

These betas are interpreted as per-SD change in phenotype per 1-SD change in 
genotype. To reflect per-SD change in phenotype per additional allele we remove
the genotype standardization (as also per SBayesRC authors 
[suggestion](https://github.com/zhilizheng/SBayesRC/issues/25#issuecomment-2138152126)):

$$
\hat\beta = \frac{Z \cdot SE}{\sqrt{2f(1-f)}}
$$

and

$$
SE = \frac{SE}{\sqrt{2f(1-f)}}
$$

The file conversions below take into account the aforementioned assumptions and
observations (indeed there are a lot of negative SEs in the data).

### Conversion of METAL to COJO

SBayesRC works with a summary statistics format supported by the parent package
GCTB called COJO:

```
METAL format
CHR   POS SNP REF ALT N   AF  SE  log10_P Z

COJO format
SNP   A1  A2  freq    b   se  p   N
```

The following R script will produce COJO versions of the METAL output:

```bash
cd $WORKSPACE/work/METAL
```

and then, write COJO format for each chromosome and genome-wide:

```r
Rscript \
  -e '{
    # Intialize output list
    cojoList <- vector("list",22)
    
    # Read files
    for (chr in seq(1,22)) {
        # METAL output filename
        mf <- paste0("metal.chr",chr,".txt")
        
        message("Reading ",mf)
        mstats <- read.delim(mf)
        
        # Intermediate variables
        # Genotype variance (for allele frequencyf=AF)
        varG <- 2 * mstats$AF * (1 - mstats$AF)
        # Denominator sqrt(N*varG)
        den <- sqrt(mstats$N * varG)
        # Standardized scale (per-SD phenotype, per-SD genotype)
        SE_sd <- 1/den
        BETA_sd <- mstats$Z * SE_sd
        
        # Assemble and write COJO file per chromosome
        cojoList[[chr]] <- data.frame(
            SNP=mstats$SNP,
            A1=mstats$ALT,
            A2=mstats$REF,
            freq=mstats$AF,
            b=BETA_sd/sqrt(varG),
            se=SE_sd/sqrt(varG),
            p=10^-mstats$log10_P,
            N=round(mstats$N)
        )
        write.table(cojoList[[chr]],file=paste0("metal_b4u_chr",chr,".ma"),
            sep="\t",quote=FALSE,row.names=FALSE)
    }
    
    # All chromosomes
    cojos <- do.call("rbind",cojoList)
    write.table(cojos,file="metal_b4u.ma",sep="\t",quote=FALSE,row.names=FALSE)
  }'
```

### Manual QC metric inspection

We check:

- SNP effect (alternative) allele frequency histograms/densities per chromosome
and in total
- Sample size per chromosome and in total
- Test various `rate2pq` thresholds against number of remaining SNPs (optional)

#### Effect allele frequency and sample size histograms/densities 

```r
library(ggplot2)

cojoList <- vector("list",22)
for (chr in seq(1,22)) {
    cj <- paste0("metal_b4u_chr",chr,".ma")
    message("Reading ",cj)
    cojoList[[chr]] <- read.delim(cj)
    cojoList[[chr]]$chr <- paste0("chr",chr)
}
cojos <- do.call("rbind",cojoList)
cojos$chr <- factor(cojos$chr,levels=paste("chr",seq(22),sep=""))

# Effect allele frequency histograms/densities per chromosome
g <- ggplot(cojos,aes(x=freq)) +
    geom_histogram(aes(y=after_stat(density)),binwidth=0.05,fill="steelblue",
        color="white") +
    geom_density(color="darkred",linewidth=1) +
    facet_wrap(~chr,ncol=4) +
    ggtitle("Allele Frequency Distribution by Chromosome") +
    xlab("\nFrequency") +
    ylab("Density\n") +
    theme_minimal() +
    theme(
        plot.title=element_text(hjust=0.5,size=16),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        strip.text.x=element_text(size=12)
    )
ggsave(filename="freq_chrom.png",plot=g,width=8,height=10)

# Effect allele frequency histograms/densities total
g <- ggplot(cojos,aes(x=freq)) +
    geom_histogram(aes(y=after_stat(density)),binwidth=0.05,fill="steelblue",
        color="white") +
    geom_density(color="darkred",linewidth=1) +
    ggtitle("Allele Frequency Distribution Genome-Wide") +
    xlab("\nFrequency") +
    ylab("Density\n") +
    theme_minimal() +
    theme(
        plot.title=element_text(hjust=0.5,size=16),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12)
    )
ggsave(filename="freq_all.png",plot=g)

# Sample size per chromosome
g <- ggplot(cojos,aes(x=N)) +
    geom_histogram(fill="steelblue",color="white") +
    facet_wrap(~chr,ncol=4) +
    ggtitle("Sample Size per Chromosome") +
    xlab("\nSize") +
    ylab("nFrequency\n") +
    theme_minimal() +
    theme(
        plot.title=element_text(hjust=0.5,size=16),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        strip.text.x=element_text(size=12)
    )
ggsave(filename="size_chrom.png",plot=g,width=8,height=10)

# Sample size total
g <- ggplot(cojos,aes(x=N)) +
    geom_histogram(fill="steelblue",color="white") +
    ggtitle("Sample Size") +
    xlab("\nSize") +
    ylab("Frequency\n") +
    theme_minimal() +
    theme(
        plot.title=element_text(hjust=0.5,size=16),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        strip.text.x=element_text(size=12)
    )
ggsave(filename="size_all.png",plot=g)
```

### Baseline PRS with R SBayesRC and 1000 genomes EUR population LD

The nature of BETTER4U per cohort summary statistics causes non-uniform and 
distant sample sizes as well as deviations from LD stability. Therefore, we run 
the QC step of SBayesRC with more relaxed threshold regarding sample size.

1. The QC and imputation steps can be done in one script as with 1000 genomes
data, not much time is required:

```bash
export OMP_NUM_THREADS=32
export RC=0.25

Rscript \
  -e '{
    library(SBayesRC)

    #WORKSPACE <- $WORKSPACE
    WORKSPACE <- Sys.getenv("WORKSPACE")
    RC <- as.numeric(Sys.getenv("RC"))

    maFile <- file.path(WORKSPACE,"work","METAL","metal_b4u.ma")
    ldFolder <- file.path(WORKSPACE,"resources","EUR_LD")
    annot <- file.path(WORKSPACE,"resources","annot_baseline2.2.txt")
    outPrefix <- file.path(WORKSPACE,"work","PRS","baseline","b4u_tgp_sbrc")
    output <- paste0(outPrefix,"_tidy.ma")

    # QC
    SBayesRC::tidy(
        mafile=maFile,
        LDdir=ldFolder,
        output=output,
        N_sd_range=4,
        log2file=TRUE
    )

    # Optional (not) imputation
    SBayesRC::impute(
        mafile=paste0(outPrefix,"_tidy.ma"),
        LDdir=ldFolder,
        output=paste0(outPrefix,"_imp.ma"),
        log2file=TRUE,
        rc=RC
    )
  }'
```

2. The tuning, MCMC iterations and effect estimations run lengthier and should
be done asynchronously, possibly with `nohup`.

```bash
export OMP_NUM_THREADS=32

nohup Rscript -e '
{
    library(SBayesRC)

    WORKSPACE <- Sys.getenv("WORKSPACE")

    ldFolder <- file.path(WORKSPACE,"resources","EUR_LD")
    annot <- file.path(WORKSPACE,"resources","annot_baseline2.2.txt")
    outPrefix <- file.path(WORKSPACE,"work","PRS","baseline","b4u_tgp_sbrc")

    # Main
    SBayesRC::sbayesrc(
        mafile=paste0(outPrefix,"_imp.ma"),
        LDdir=ldFolder,
        outPrefix=paste0(outPrefix,"_prs"),
        annot=annot,
        log2file=TRUE
    )
}
' > $WORKSPACE/work/PRS/baseline/b4u_tgp_sbrc_prs.log 2>&1 &
```

### Baseline PRS with R SBayesRC and built-in LD

1. We run with more relaxed thresholds due to BETTER4U data nature. Also, 
imputation takes longer, so we run in the background with `nohup`.

```bash
export OMP_NUM_THREADS=32
export RC=0.5

nohup Rscript -e '
{
    library(SBayesRC)

    WORKSPACE <- Sys.getenv("WORKSPACE")
    RC <- as.numeric(Sys.getenv("RC"))

    maFile <- file.path(WORKSPACE,"work","METAL","metal_b4u.ma")
    ldFolder <- file.path(WORKSPACE,"resources","ukbEUR_Imputed")
    annot <- file.path(WORKSPACE,"resources","annot_baseline2.2.txt")
    outPrefix <- file.path(WORKSPACE,"work","PRS","baseline","b4u_ukb_sbrc")
    output <- paste0(outPrefix,"_tidy.ma")

    # QC
    SBayesRC::tidy(
        mafile=maFile,
        LDdir=ldFolder,
        output=output,
        N_sd_range=4,
        log2file=TRUE
    )

    # Optional (not) imputation
    SBayesRC::impute(
        mafile=paste0(outPrefix,"_tidy.ma"),
        LDdir=ldFolder,
        output=paste0(outPrefix,"_imp.ma"),
        log2file=TRUE,
        rc=RC
    )
}
' > $WORKSPACE/work/PRS/baseline/b4u_ukb_sbrc_imp.log 2>&1 &
```

2. The tuning, MCMC iterations and effect estimations run lengthier and should
be done asynchronously, possibly with `nohup`.

```bash
export OMP_NUM_THREADS=32

nohup Rscript -e '
{
    library(SBayesRC)

    WORKSPACE <- Sys.getenv("WORKSPACE")

    ldFolder <- file.path(WORKSPACE,"resources","ukbEUR_Imputed")
    annot <- file.path(WORKSPACE,"resources","annot_baseline2.2.txt")
    outPrefix <- file.path(WORKSPACE,"work","PRS","baseline","b4u_ukb_sbrc")

    # Main
    SBayesRC::sbayesrc(
        mafile=paste0(outPrefix,"_imp.ma"),
        LDdir=ldFolder,
        outPrefix=paste0(outPrefix,"_prs"),
        annot=annot,
        log2file=TRUE
    )
}
' > $WORKSPACE/work/PRS/baseline/b4u_ukb_sbrc_prs.log 2>&1 &
```
2814047
2814105
### Baseline PRS with GCTB and 1000 genomes EUR population LD

The QC step in GCTB SBayesRC is fixed and uses only _N_.

1. The QC and imputation steps can be done in one script as with 1000 genomes
data, not much time is required:


```bash
THREADS=32

gctb \
  --ldm-eigen $WORKSPACE/resources/EUR_LD \
  --gwas-summary $WORKSPACE/work/METAL/metal_b4u.ma \
  --impute-summary \
  --out $WORKSPACE/work/PRS/baseline/b4u_tgp_gctb \
  --thread $THREADS
```

2. With 1000 genomes EUR LD, some SNPs are missing from the functional 
annotations and if we don't fix, error occurs.

```bash
cut -f1 $WORKSPACE/resources/annot_baseline2.2.txt \
   > $WORKSPACE/resources/annot_snps.txt

Rscript \
  -e '{
    WORKSPACE <- Sys.getenv("WORKSPACE")

    ref <- read.delim(file.path(WORKSPACE,"resources","annot_snps.txt"))
    gwa <- read.table(file.path(WORKSPACE,"work","PRS","baseline",
        "b4u_tgp_gctb.imputed.ma"),header=TRUE)
    reqs <- setdiff(gwa$SNP,ref$SNP)
    tmp <- read.delim(file.path(WORKSPACE,"resources","annot_baseline2.2.txt"),
        nrow=5)
    out <- matrix(0,nrow=length(reqs),ncol=ncol(tmp)-1)
    colnames(out) <- names(tmp)[2:ncol(tmp)]
    out[,1] <- 1
    out <- data.frame(reqs,out)
    names(out)[1] <- "SNP"
    write.table(out,file=file.path(WORKSPACE,"resources","annot_add.txt"),
        quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
  }'

cat $WORKSPACE/resources/annot_baseline2.2.txt $WORKSPACE/resources/annot_add.txt \
 > $WORKSPACE/resources/annot_expanded.txt
```

3. GCTB SBayesRC (tuning, MCMC iterations and effect estimations) is now ready 
to run (<1h)

```bash
gctb \
  --ldm-eigen $WORKSPACE/resources/EUR_LD \
  --gwas-summary $WORKSPACE/work/PRS/baseline/b4u_tgp_gctb.imputed.ma \
  --sbayes RC \
  --annot $WORKSPACE/resources/annot_expanded.txt \
  --out $WORKSPACE/work/PRS/baseline/b4u_tgp_gctb \
  --thread $THREADS
```

### Baseline PRS with GCTB and built-in LD

The QC step in GCTB SBayesRC is fixed and uses only _N_.

1. The QC and imputation steps can be done in one script as with 1000 genomes
data, not much time is required:

```bash
THREADS=32

gctb \
  --ldm-eigen $WORKSPACE/resources/ukbEUR_Imputed \
  --gwas-summary $WORKSPACE/work/METAL/metal_b4u.ma \
  --impute-summary \
  --out $WORKSPACE/work/PRS/baseline/b4u_ukb_gctb \
  --thread $THREADS
```

2. GCTB SBayesRC (tuning, MCMC iterations and effect estimations) is now ready 
to run (requires ~6-7h, better with `nohup`)

```bash
nohup gctb \
  --ldm-eigen $WORKSPACE/resources/ukbEUR_Imputed \
  --gwas-summary $WORKSPACE/work/PRS/baseline/b4u_ukb_gctb.imputed.ma \
  --sbayes RC \
  --annot $WORKSPACE/resources/annot_baseline2.2.txt \
  --out $WORKSPACE/work/PRS/baseline/b4u_ukb_gctb \
  --thread $THREADS \
  > $WORKSPACE/work/PRS/baseline/b4u_ukb_gctb.log 2>&1 &
```

## Baseline PRS derivation with PRS-CS

### Conversion of METAL to PRS-CS expected format

PRS-CS(x) works with the following summary statistics format:

```
METAL format
CHR   POS SNP REF ALT N   AF  SE  log10_P Z

PRS-CSx format
SNP          A1   A2   BETA      SE
```

The following R script will produce PRS-CSx versions of the METAL output. The
same assumptions regarding METAL output and betas, SEs as in the case of 
SBayesRC COJO files apply:

```bash
cd $WORKSPACE/work/METAL
```

and then, write PRS-CS(x) format for each chromosome and genome-wide:

```r
Rscript \
  -e '{
    # Intialize output list
    prscxList <- vector("list",22)
    
    # Read files and convert
    for (chr in seq(1,22)) {
        # METAL output filename
        mf <- paste0("metal.chr",chr,".txt")
        
        message("Reading ",mf)
        mstats <- read.delim(mf)
        
        # Intermediate variables
        # Genotype variance (for allele frequencyf=AF)
        varG <- 2 * mstats$AF * (1 - mstats$AF)
        # Denominator sqrt(N*varG)
        den <- sqrt(mstats$N * varG)
        # Standardized scale (per-SD phenotype, per-SD genotype)
        SE_sd <- 1/den
        BETA_sd <- mstats$Z * SE_sd
        
        # Assemble and write PRS-CS files per chromosome
        prscxList[[chr]] <- data.frame(
            SNP=mstats$SNP,
            A1=mstats$ALT,
            A2=mstats$REF,
            BETA=BETA_sd/sqrt(varG),
            SE=SE_sd/sqrt(varG),
            p=10^-mstats$log10_P,
            N=round(mstats$N)
        )
        write.table(prscxList[[chr]],file=paste0("metal_b4u_chr",chr,".csx"),
            sep="\t",quote=FALSE,row.names=FALSE)
    }
    
    # All chromosomes
    prscxs <- do.call("rbind",prscxList)
    write.table(prscxs,file="metal_b4u.csx",sep="\t",quote=FALSE,
        row.names=FALSE)
  }'
```

### Baseline PRS with PRS-CS and 1000 genomes EUR population LD

PRS-CS requires a PLINK BIM file for the variants in the validation populatyion.
BED files are not required as PRS-CS can operate on summary statistics. In this
case the BIM file is used for overlap with the 1000 genomes LD. We make sure 
that the BIM file lives on its own:

```bash
mkdir -p $WORKSPACE/work/HUABB/only_bim && cd $WORKSPACE/work/HUABB/only_bim
cp ../HUA_unrelated_dbsnp ./
cd ..
```

PRS-CS runs in a single step. It may take some time to complete, so it's better
to use `nohup`:

```bash
# --n_gwas can be the median METAL sample size
nohup python $WORKSPACE/resources/PRScs/PRScs.py \
  --ref_dir=$WORKSPACE/resources/PRScsxLD/ldblk_1kg_eur \
  --bim_prefix=$WORKSPACE/work/HUABB/only_bim/HUA_unrelated_dbsnp \
  --sst_file=$WORKSPACE/work/METAL/metal_b4u.csx \
  --n_gwas=232593 \
  --out_dir=$WORKSPACE/work/PRS/baseline/b4u_tgp_prscs \
  > $WORKSPACE/work/PRS/baseline/b4u_tgp_prscs.log 2>&1 &
# 2817564
```

Then simply concatenate the per chromosome results:

```bash
cd $WORKSPACE/work/PRS/baseline
for i in {1..22}; 
do
  cat b4u_tgp_prscs_pst_eff_a1_b0.5_phiauto_chr${i}.txt
done > b4u_tgp_prscs_pst_eff_a1_b0.5_phiauto.txt
```

## Baseline PRS calculations and evaluation

### Basic filtering of R SBayesRC derived PRS

The R SBayesRC PRS file has the format:

```
SNP A1  BETA    SE  PIP BETAlast
rs12565286  C   1.1512915188263e-05 0.00114985060792772 0.00349998474121094 0
rs3115850   C   7.64806479735764e-06    0.00053908843748555 0.00150001049041748 0
rs2286139   T   0   0   0   0
...
```

We keep the SNP, A1, BETA, SE (for later introducing noise), PIP and remove SNPs
with zero effect. We process both 1000 genomes EUR LD and built-in UKB derived 
PRS:

```
# 1000 genomes
awk 'NR==1 || $3 != 0 { print $1"\t"$2"\t"$3"\t"$4"\t"$5 }' \
  $WORKSPACE/work/PRS/baseline/b4u_tgp_sbrc_prs.txt > \
  $WORKSPACE/work/PRS/baseline/b4u_tgp_sbrc_prs.prs

# Built-in
awk 'NR==1 || $3 != 0 { print $1"\t"$2"\t"$3"\t"$4"\t"$5 }' \
  $WORKSPACE/work/PRS/baseline/b4u_ukb_sbrc_prs.txt > \
  $WORKSPACE/work/PRS/baseline/b4u_ukb_sbrc_prs.prs
```

### Basic filtering of GCTB SBayesRC derived PRS

The GCTB SBayesRC PRS file has the format:

```
 Index                 Name  Chrom     Position     A1     A2        A1Frq     A1Effect           SE VarExplained          PEP          Pi1          Pi2          Pi3          Pi4          Pi5            PIP
     1           rs12565286      1       721290      C      G     0.045726     0.000000     0.000000 0.000000e+00     0.000000     0.994805     0.004617     0.000549     0.000029     0.000000      0.0051949
     2            rs3094315      1       752566      A      G     0.839960     0.000140     0.002029 4.040904e-08     0.005000     0.995331     0.004260     0.000391     0.000017     0.000000     0.00466907
...
```

We keep the Name, A1, A1Effect, SE (for later introducing noise), PIP and 
remove SNPs with zero effect. We process both 1000 genomes EUR LD and built-in 
UKB derived PRS:

```
# 1000 genomes
awk 'NR==1 { print "SNP\tA1\tBETA\tSE\tPIP" } 
     NR>1 && $8 != 0 { print $2"\t"$5"\t"$8"\t"$9"\t"$17 }' \
  $WORKSPACE/work/PRS/baseline/b4u_tgp_gctb.snpRes > \
  $WORKSPACE/work/PRS/baseline/b4u_tgp_gctb.snpRes.prs

# Built-in
awk 'NR==1 { print "SNP\tA1\tBETA\tSE\tPIP" } 
     NR>1 && $8 != 0 { print $2"\t"$5"\t"$8"\t"$9"\t"$17 }' \
  $WORKSPACE/work/PRS/baseline/b4u_ukb_gctb.snpRes > \
  $WORKSPACE/work/PRS/baseline/b4u_ukb_gctb.snpRes.prs
```

### Baseline PRS evaluation for both SBayesRC and PRS-CS

We developed a small R library to sanitize (allele flipping etc.) and evaluate 
the derived PRS based on simple metrics such as:

- R<sup>2</sup> of a simple regression model fit with the phenotype (BMI etc) as 
a function of the covariates (sex, age, etc) excluding the PRS as a covariate
(the _null_ model)
- R<sup>2</sup> of a simple regression model fit with the phenotype (BMI etc) as 
a function of the covariates (sex, age, etc) including the PRS as a covariate
(the _full_ model)
- The difference between R<sup>2</sup> of the _full_ model and the _null_ model
(PRS R<sup>2</sup>)
- The statistical significance of the inclusion of the PRS in the model (PRS
p-value)
- Correlations between the predicted phenotype values and the actual phenotype
values when using the _null_ and _full_ models
- Correlation between the PRS and the phenotype
- Statistical significance of the correlation between the PRS and the phenotype

Below, we evaluate only based on PRS R<sup>2</sup> but other values can be
reported.

Furthermore, SBayesRC results can be subset and filtered based on i) betas, ii)
posterior SNP probabilities (PIP) and PRS-CS results can be filtered based on
betas. In order to keep the number of SNPs in the PRS as small as possible (so
that we obtain good coverage across cohorts) we implemented a simple grid search
algorithm which evaluates the chosen metric for several beta and PIP (where
available) values. The functions accept more arguments than those mentioned
below, the code is available and can be explored. PLINK is required to be 
present in the system.

The whole process is outlined in the R chunks below:

```r
source("evalfuns.R")

WORKSPACE <- Sys.getenv("WORKSPACE")

# Genetic data file (HUABB)
genoBase <- file.path(WORKSPACE,"work","HUABB","HUA_unrelated_dbsnp")
# Covariates file (sex, age, age^2, PCs) and trait
covFile <- file.path(WORKSPACE,"work","HUABB","huabb","HUA_covariates.txt")
trait <- "bmi"

# R SBayesRC with 1000 genomes LD panel
prsFile <- file.path(WORKSPACE,"work","PRS","baseline","b4u_tgp_sbrc_prs.prs")

# Grid search flow
# 1. Sanitize PRS
sanFile <- sanitizePrs(prsFile,genoBase)

# 2. Grid search (absolute values of beta)
obj <- gridSearch(prsFile=sanFile,covFile=covFile,trait=trait,genoBase=genoBase,
    rc=0.2)
#Best PRS according to prs_r2 found at:
#  abs(BETA) > 1e-05
#  PIP > 0
# i = 7, j = 1

# 3. Plot results
p <- gridSearchPlot(obj$metrics,"prs_r2",i=7,j=1)


# R SBayesRC with built-in LD panel
prsFile <- file.path(WORKSPACE,"work","PRS","baseline","b4u_ukb_sbrc_prs.prs")

# 1. Sanitize PRS
sanFile <- sanitizePrs(prsFile,genoBase)

# 2. Grid search (absolute values of beta)
obj <- gridSearch(prsFile=sanFile,covFile=covFile,trait=trait,genoBase=genoBase,
    rc=0.2)
#Best PRS according to prs_r2 found at:
#  abs(BETA) > 1e-06
#  PIP > 0
# i = 6, j = 1

# 3. Plot results
p <- gridSearchPlot(obj$metrics,"prs_r2",i=6,j=1)


# GCTB SBayesRC with 1000 genomes LD panel
prsFile <- file.path(WORKSPACE,"work","PRS","baseline","b4u_tgp_gctb.snpRes.prs")

# 1. Sanitize PRS
sanFile <- sanitizePrs(prsFile,genoBase)

# 2. Grid search (absolute values of beta)
obj <- gridSearch(prsFile=sanFile,covFile=covFile,trait=trait,genoBase=genoBase,
    rc=0.2)
#Best PRS according to prs_r2 found at:
#  abs(BETA) > 1e-04
#  PIP > 0.1
# i = 8, j = 6

# 3. Plot results
p <- gridSearchPlot(obj$metrics,"prs_r2",i=8,j=6)


# GCTB SBayesRC with built-in LD panel
prsFile <- file.path(WORKSPACE,"work","PRS","baseline","b4u_ukb_gctb.snpRes.prs")

# 1. Sanitize PRS
sanFile <- sanitizePrs(prsFile,genoBase)

# 2. Grid search (absolute values of beta)
obj <- gridSearch(prsFile=sanFile,covFile=covFile,trait=trait,genoBase=genoBase,
    rc=0.2)
#Best PRS according to prs_r2 found at:
#  abs(BETA) > 1e-05
#  PIP > 0.01
# i = 8, j = 6

# 3. Plot results
p <- gridSearchPlot(obj$metrics,"prs_r2",i=8,j=6)


# PRS-CS with built-in EUR 1000 genomes panel
prsFile <- file.path(WORKSPACE,"work","PRS","baseline","b4u_tgp_prscs_pst_eff_a1_b0.5_phiauto.txt")

# 1. Sanitize PRS
sanFile <- sanitizePrs(prsFile,genoBase,from="prscs")

# 2. Grid search (absolute values of beta)
obj <- gridSearch(prsFile=sanFile,covFile=covFile,trait=trait,genoBase=genoBase,
    rc=0.2)
#Best PRS according to prs_r2 found at:
#  abs(BETA) > 1e-06
#  PIP > 0
# i = 6, j = 1
```

Based on the aforementioned evaluations, it seems that PRS-CS is the choice to
go for BMI (and subsequently weight/BMI change).

## PRS bootstrapping with PRS-CS

### Theoretical framework

We have SNP effects from meta-analyzed summary statistics from all the BETTER4U 
cohosrts. We seek to add Gaussian noise to the beta effects such that the 
resulting betas fall within their 95% estimated confidence interval.

We have:

* $\hat{\beta}_j$ : estimated effect size for SNP<sub>j</sub>
* $SE_j$ : standard error of $\hat{\beta}_j$

Then the 95% confidence interval for $\hat{\beta}_j$ is:

$$
\left[ \hat{\beta}_j - 1.96 \cdot SE_j,\ \hat{\beta}_j + 1.96 \cdot SE_j \right]
$$

We add Gaussian noise to simulate perturbation:

$$
\tilde{\beta}_j = \hat{\beta}_j + \epsilon_j,\quad \text{where } \epsilon_j \sim \mathcal{N}(0, \lambda \cdot SE_j^2)
$$

Now, we want to estimate the probability the perturbed beta stays in the 95% CI.
We want:

$$
P\left( \tilde{\beta}_j \in [\hat{\beta}_j \pm 1.96 \cdot SE_j] \right) = P\left( |\epsilon_j| < 1.96 \cdot SE_j \right) = 0.95
$$

With error rescaling:

$$
Z_j = \frac{\epsilon_j}{SE_j} \sim \mathcal{N}(0, \lambda)
$$

the inequality becomes:

$$
P(|Z_j| < 1.96) = 0.95
$$

We are now solving for λ and looking for the value of λ such that:

$$
P(|Z| < 1.96) = 0.95 \quad \text{where } Z \sim \mathcal{N}(0, \lambda)
$$

If Z were from a **standard normal** $\mathcal{N}(0,1)$, this would be true by 
definition. But here the **variance** is λ, so the tails are fatter (for λ > 1) 
or thinner (for λ < 1).

So we solve:

$$
\Phi\left( \frac{1.96}{\sqrt{\lambda}} \right) - \Phi\left( -\frac{1.96}{\sqrt{\lambda}} \right) = 0.95
$$

This is the cumulative probability from −1.96 to +1.96 under a normal 
distribution with variance λ. If we solve the equation numerically, we get:

$$
\lambda \approx 1.0
$$

This makes intuitive sense because if the perturbation variance equals the 
original SE<sup>2</sup> (i.e. λ = 1), then:

$$
\epsilon_j \sim \mathcal{N}(0, SE_j^2)
\Rightarrow Z_j = \epsilon_j / SE_j \sim \mathcal{N}(0, 1)
\Rightarrow P(|Z_j| < 1.96) \approx 0.95
$$

So, **λ = 1** means we are adding noise of the same scale as the uncertainty in 
the original beta estimates.

### Generation of perturbed summary statistics

The file `$WORKSPACE/work/METAL/metal_b4u.csx` contains the meta-analysis
summary statistics in PRS-CS format. We will use betas and SEs in this file to
generate pertrurbed summary statistics. The following script will introduce 
noise in betas within the acceptable CIs. The process will be repeated 1000
times to generate 1000 perturbed datasets.

```
source("evalfuns.R")

WORKSPACE <- Sys.getenv("WORKSPACE")

# Read in summary statistics
gwas <- read.delim(file.path(WORKSPACE,"work","METAL","metal_b4u.csx"))


```


The following script will execute PRS-CS 1000 times

PRS-CS runs in a single step. It may take some time to complete, so it's better
to use `nohup`:

```bash
# --n_gwas can be the median METAL sample size
nohup python $WORKSPACE/resources/PRScs/PRScs.py \
  --ref_dir=$WORKSPACE/resources/PRScsxLD/ldblk_1kg_eur \
  --bim_prefix=$WORKSPACE/work/HUABB/only_bim/HUA_unrelated_dbsnp \
  --sst_file=$WORKSPACE/work/METAL/metal_b4u.csx \
  --n_gwas=232593 \
  --out_dir=$WORKSPACE/work/PRS/baseline/b4u_tgp_prscs \
  > $WORKSPACE/work/PRS/baseline/b4u_tgp_prscs.log 2>&1 &
# 2614263
```

Then simply concatenate the per chromosome results:

```bash
cd $WORKSPACE/work/PRS/baseline
for i in {1..22}; 
do
  cat b4u_tgp_prscs_pst_eff_a1_b0.5_phiauto_chr${i}.txt
done > b4u_tgp_prscs_pst_eff_a1_b0.5_phiauto.txt
```


#### Notes

1. Manual calculation of the `rate2pq` metric

```r
D <- 2 * z$freq * (1 - z$freq) * z$N
varps <- D*(z$N * z$se^2 + z$b^2)/z$N
vary <- median(varps)
indic <- sqrt( 2 * z$freq * (1 - z$freq) / vary * (z$N*z$se^2 + z$b^2))
```

2. Transformations to original beta scale (not really helping)

```r
b = mstats$Z*mstats$SE / sqrt(2*mstats$AF*(1 - mstats$AF))
se = mstats$SE / sqrt( 2 * mstats$AF * (1 - mstats$AF))
```

3. LD with 1000 genomes gives slightly more SNPs after QC

4. Imputation with built-in LD

Place the following in a temporary R script in the same folder as the data
or somewhere else as long as the WORKSPACE directory is defined.

```bash
cat <<EOF > b4uukb_impute.R
library(SBayesRC)

WORKSPACE <- "/media/storage3/playground/b4uprs"

maFile <- file.path(WORKSPACE,"/work/METAL/metal_b4u.ma")
ldFolder <- file.path(WORKSPACE,"resources/ukbEUR_Imputed")
annot <- file.path(WORKSPACE,"resources/annot_baseline2.2.txt")
outPrefix <- file.path(WORKSPACE,"work/METAL/b4uukb")

SBayesRC::impute(
    mafile=paste0(outPrefix,"_tidy.ma"),
    LDdir=ldFolder,
    output=paste0(outPrefix,"_imp.ma")
)
EOF

nohup Rscript b4uukb_impute.R > b4uukb_impute.log &
```

5. We have 4 approaches:
 - Baseline PRS with R SBayesRC and 1000 genomes LD
 - Baseline PRS with R SBayesRC and built-in LD
 - Baseline PRS with GCTB SBayesRC and 1000 genomes LD
 - Baseline PRS with GCTB SBayesRC and built-in LD
We will decide on what to iterate and perturb based on the best R2 with HUA toy 
or biobank. The main difference is the QC step

6. There are too many SNPs in the output, especially with the built-in LD. We
need to test with a grid of `BETA` and `PIP` thresholds to evaluate R<sup>2</sup>
with the BMI phenotype in HUABB. The coverage should also be an output from the
gridding.

7. Some facts for further steps

The baseline PRS derivation process leads to 4 baseline PRS. We now need to 
apply them in order to decide which one will go through bootstrapping.

Two output formats:
- SBayesRC
- GCTB

**SBayesRC**

```
SNP A1  BETA    SE  PIP BETAlast
rs12565286  C   1.1512915188263e-05 0.00114985060792772 0.00349998474121094 0
rs3115850   C   7.64806479735764e-06    0.00053908843748555 0.00150001049041748 0
rs2286139   T   0   0   0   0
```

**GCTB**

```
Index                 Name  Chrom     Position     A1     A2        A1Frq     A1Effect           SE VarExplained          PEP          Pi1          Pi2          Pi3          Pi4          Pi5            PIP
1           rs12565286      1       721290      C      G     0.045726    -0.000016     0.000224 1.605210e-10     0.000000     0.995424     0.003256     0.001219     0.000102     0.000000     0.00457633
```

These need to be preprocessed and scored against PLINK files.

Regarding SBayesRC output:
- Remove SNPs with 0 effect
- Keep SNP, A1, BETA, SE (for introducing noise)

Regarding GCTB output
- Remove SNPs with 0 effect
- Keep Name, Chrom, Position, A1, A2, A1Effect, SE (for introducing noise)

Regarding PLINK files (handled in small R library developed):
- Conversion from dbSNP annotated VCF while keeping allele order
- Ensure allele matching/flipping
- Do not switch in PLINK, rather change betas in score files

If we have VCF we need to convert to PLINK for scoring. In the case of one VCF
per chromosome:

```
for CHR in `seq 1 22`
do
  plink \
    --vcf COHORT_${CHR}.vcf.gz \
    --keep-allele-order \
    --vcf-idspace-to _ \
    --const-fid \
    --allow-extra-chr 0 \
    --split-x b37 no-fail \
    --make-bed \
    --out COHORT &
done
```

or in the case of one merged VCF:

```
plink \
  --vcf COHORT.vcf.gz \
  --keep-allele-order \
  --vcf-idspace-to _ \
  --const-fid \
  --allow-extra-chr 0 \
  --split-x b37 no-fail \
  --make-bed \
  --out COHORT
```

8. Preprocessing of HUABB

```bash
bcftools query -l HUA_unrelated_dbsnp.vcf.gz | cut -d'_' -f1 > proper_names.txt

bcftools annotate -x ID HUA_unrelated_dbsnp.vcf.gz -Oz \
  -o HUA_unrelated_tmp.vcf.gz

bcftools reheader -s proper_names.txt \
  -o HUA_unrelated.vcf.gz  HUA_unrelated_tmp.vcf.gz

rm HUA_unrelated_tmp.vcf.gz

nohup bash -c '
java -jar ../../resources/snpEff/SnpSift.jar annotate -noInfo -noLog \
    -id ../../resources/dbsnp/dbSNP155.vcf.gz \
    HUA_unrelated.vcf.gz |
  bcftools view --include 'ID != "."' \
    --output-type z \
    --output HUA_unrelated_dbsnp.vcf.gz' > annotate.log 2>&1 &

bcftools view -H HUA_unrelated_dbsnp_pre.vcf.gz | cut -f1,2 | sort | uniq -d | \
  bcftools view -T ^<(cat) -Oz -o HUA_unrelated_dbsnp.vcf.gz HUA_unrelated_dbsnp_pre.vcf.gz

plink \
  --vcf HUA_unrelated_dbsnp.vcf.gz \
  --keep-allele-order \
  --vcf-idspace-to _ \
  --const-fid \
  --allow-extra-chr 0 \
  --split-x b37 no-fail \
  --make-bed \
  --out HUA_unrelated_dbsnp
```

9. Potential PRS-CS(x) call

```bash
python PRScsx.py \
  --ref_dir=$WORKSPACE/resources/PRScsxLD \
  --bim_prefix=$WORKSPACE/work/HUABB/only_bim/HUA_unrelated_dbsnp \
  --sst_file=$WORKSPACE/work/METAL/metal_b4u.csx \
  --n_gwas=230146 \
  --pop=EUR \
  --out_dir=$WORKSPACE/work/PRS/baseline/b4u_tgp_prscs \
  --out_name=b4u_tgp_prscs
```

```bash
nohup python $WORKSPACE/resources/PRScs/PRScs.py \
 --ref_dir=$WORKSPACE/resources/PRScsxLD/ldblk_1kg_eur \
 --bim_prefix=$WORKSPACE/work/HUABB/only_bim/HUA_unrelated_dbsnp \
 --sst_file=$WORKSPACE/work/METAL/metal_b4u.csx \
 --n_gwas=230146 \
 --out_dir=b4u_tgp_prscs &
```

10. Initial COJO creation prior to re-examination of METAL output

```r
Rscript \
  -e '{
    cojoList <- vector("list",22)
    for (chr in seq(1,22)) {
        mf <- paste0("metal.chr",chr,".txt")
        message("Reading ",mf)
        mstats <- read.delim(mf)
        cojoList[[chr]] <- data.frame(
            SNP=mstats$SNP,
            #A1=mstats$REF,
            #A2=mstats$ALT,
            A1=mstats$ALT,
            A2=mstats$REF,
            freq=mstats$AF,
            b=mstats$Z*mstats$SE,
            #b=mstats$Z*mstats$SE / sqrt(2*mstats$AF*(1 - mstats$AF)),
            se=mstats$SE,
            #se=mstats$SE / sqrt( 2 * mstats$AF * (1 - mstats$AF)),
            p=10^-mstats$log10_P,
            N=round(mstats$N)
        )
        write.table(cojoList[[chr]],file=paste0("metal_b4u_chr",chr,".ma"),
            sep="\t",quote=FALSE,row.names=FALSE)
    }
    cojos <- do.call("rbind",cojoList)
    write.table(cojos,file="metal_b4u.ma",sep="\t",quote=FALSE,row.names=FALSE)
  }'
```

11. Initual SBayesRC running with wrong betas and relaxed thresholds

```bash
export OMP_NUM_THREADS=32
export RC=0.25

Rscript \
  -e '{
    library(SBayesRC)

    #WORKSPACE <- $WORKSPACE
    WORKSPACE <- Sys.getenv("WORKSPACE")
    RC <- as.numeric(Sys.getenv("RC"))

    maFile <- file.path(WORKSPACE,"work","METAL","metal_b4u.ma")
    ldFolder <- file.path(WORKSPACE,"resources","EUR_LD")
    annot <- file.path(WORKSPACE,"resources","annot_baseline2.2.txt")
    outPrefix <- file.path(WORKSPACE,"work","PRS","baseline","b4u_tgp_sbrc")
    output <- paste0(outPrefix,"_tidy.ma")

    # QC
    SBayesRC::tidy(
        mafile=maFile,
        LDdir=ldFolder,
        output=output,
        freq_thresh=0.3,
        N_sd_range=4,
        rate2pq=0.8,
        log2file=TRUE
    )

    # Optional (not) imputation
    SBayesRC::impute(
        mafile=paste0(outPrefix,"_tidy.ma"),
        LDdir=ldFolder,
        output=paste0(outPrefix,"_imp.ma"),
        log2file=TRUE,
        rc=RC
    )
  }'
```
