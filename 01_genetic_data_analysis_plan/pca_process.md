Calculation of main Principal Components for BETTER4U
================================================================================

### Authors

* Panagiotis Moulos (moulos@fleming.gr)
* Anders Eriksson (anders.eriksson@ut.ee)

The following part for performing PCA on 1000 genomes data is based on
[this](https://www.biostars.org/p/335605/) tutorial from 
[Kevin Blighe](https://github.com/kevinblighe) in Biostars.

# Prerequisites

* VCF files per chromosomes from the latest 1000 genomes project website
* [PLINK 1.90](https://www.cog-genomics.org/plink/)
* [bcftools](https://samtools.github.io/bcftools/bcftools.html)
* [FlashPCA2](https://github.com/gabraham/flashpca)

# Process

## 1. Download the files as VCF.gz (and tab-indices)


```
PREFIX="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.CHR" ;
SUFFIX=".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" ;

for CHR in `seq 1 22`; 
do
    #echo ${PREFIX}${CHR}${SUFFIX}
    #echo ${PREFIX}${CHR}${SUFFIX}".tbi"
    wget ${PREFIX}${CHR}${SUFFIX} ${PREFIX}${CHR}${SUFFIX}".tbi";
done
```

## 2. Download the GRCh37/hg19 reference genome

For better data integority, we download the human reference genome version used
by the consortium:

```
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz ;
```

## 3. Convert the 1000 Genomes files to BCF

In the following we are creating a first set of BCF files from 1000 genomes data
which is suitable for further dowstreams selection of variants for PCA for the
BETTER4U cohorts:

* Include only bi-allelic SNPs with MAF > 0.05 (1st pipe)
* Ensure that multi-allelic calls are split and that indels are left-aligned 
compared to reference genome (2nd pipe)
* Set the VCF ID field to the value: `CHROM:POS` (3rd pipe, this is suitable for
using later with HRC variants as only unique SNPs are contained)
* Remove duplicates (4th pipe)

Notes:

- `-I +'%CHROM:%POS'` means that unset IDs will be set to 
`CHROM:POS*`
- `-x ID -I +'%CHROM:%POS:%REF:%ALT'` first erases the current ID and then sets 
it to `CHROM:POS*`

```
for CHR in `seq 1 22`; 
do
  echo "Converting $CHR"
  
  bcftools view --include 'INFO/AF > 0.05' --types snps \
    --min-alleles 2 --max-alleles 2 \
    ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz | \
  bcftools norm --multiallelics -any --check-ref w \
    --fasta-ref human_g1k_v37.fasta | \
  bcftools annotate --remove ID --set-id +'%CHROM:%POS' | \
  bcftools norm --output-type b --rm-dup both \
    --output ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf;
done
```

Then index:

```
for CHR in `seq 1 22`; 
do
  echo "Indexing $CHR"
  
  bcftools index \    
    ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf;
done
```

And get a list of the remaining variants in `CHROM:POS` format to later 
intersect with the HRC panel variants:

```
for CHR in `seq 1 22`; 
do
  echo "Getting ids for $CHR"
  
  bcftools query --format '%CHROM:%POS' \    
    ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf \
    > ALL.chr${CHR}.variant.ids;
done
```

The above command will probably take some time, so they can be put in a shell
script and executed with `nohup`. This script should look like the following 
(named `1000g2bcf.sh`):

```
#!/bin/bash

for CHR in `seq 1 22`;
do
  echo "Converting $CHR"

  bcftools view --include 'INFO/AF > 0.05' --types snps \
    --min-alleles 2 --max-alleles 2 \
    ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz | \
  bcftools norm --multiallelics -any --check-ref w \
    --fasta-ref human_g1k_v37.fasta | \
  bcftools annotate --remove ID --set-id +'%CHROM:%POS' | \
  bcftools norm --output-type b --rm-dup both \
    --output ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf &
done

wait

for CHR in `seq 1 22`;
do
  echo "Indexing $CHR"

  bcftools \
    index ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf &
done

wait

for CHR in `seq 1 22`;
do
  echo "Getting ids for $CHR"

  bcftools query --format '%CHROM:%POS' \
    ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf > \
      ALL.chr${CHR}.variant.ids &
done
```

and the command to execute it:

```
nohup sh 1000g2bcf.sh > 1000g2bcf.log &
```

## 4. Prepare a list of common SNPs accross partners

This steps compiles a list of good quality imputed SNPs from all the partner
cohorts based on information gathered from each partner for the HRC imputed
SNPs using a script provided in Appendix I. The creation of this list is 
outlined in the following.

### List of partners with genetic data:

The partners that will provide genetic data are:

- HUA
- REGIONH
- MUW
- UCY
- UTARTU
- HMGU
- VUA
- TAUH
- QIMR
- BIB
- WIS

Most partners provided joint array data (possibly different platforms) per 
study. Some provided imputed data per study and imputed platform. Every imputed
dataset, either joint or per platform, comprising imputed chromosomes 1 to 22
was treated as a separate set per partner. Based on this assumption, the 
following table summarizes what we received:


| Partner | Cohort             | Datasets |
| ------- | ------------------ | -------- |
| UTARTU  | EstBB              | 1        |
| HMGU    | GINIplus&LISA      | 1        |
| TAUH    | ANGES              | 1        |
| TAUH    | FINCANVAS_HCE      | 1        |
| TAUH    | FINCANVAS_MC       | 1        |
| TAUH    | YFS                | 1        |
| REGIONH | COPSAC2000         | 1        |
| REGIONH | COPSAC2010         | 1        |
| REGIONH | COPSAC2010_PARENTS | 1        |
| MUW     | HPFS               | 6        |
| MUW     | NHS                | 5        |
| MUW     | NHS2               | 4        |
| MUW     | PHS                | 3        |
| HUA     | MANOLIS            | 1        |
| HUA     | POMAC              | 1        |
| HUA     | NAFLD              | 1        |
| HUA     | TEENAGE            | 1        |
| HUA     | THISEAS            | 1        |
| HUA     | IMPROVE            | 1        |
| HUA     | OSTEOS             | 1        |
| VUA     | NTR                | 1        |


### Common SNP extraction process

The data were retrieved and one directory for each partner was created. Next,
data for each cohort and dataset were placed inside each partner directory
without making distinctions (for the purpose of this work) between a cohort with
one dataset and a cohort with multiple datasets. For example, UTARTU provided 
one cohort consisting of 1 (merged) dataset, TAUH provided four cohorts 
consisting of one dataset each, while MUW provided four cohorts comprising 
multiple dataset each. Each dataset was placed in the directory of each partner
without any further processing (e.g. merging datasets from MUW), therefore, the
directory of UTARTU contained one subdirectory, the directory of TAUH contained
four subdirectories while the directory of MUW contained 18 subdirectories. 

Based on the above assumptions, the following R script:

* Traverses each partner directory and reads the subdirectory contents (partner
datasets)
* Filters each dataset according to the Minimac INFO or Beagle R2 score 
(cutoff: 0.3)
* Gathers all the good quality SNPs from all individual partner datasets
* Find the common good SNPs from all partners for each chromosome and writes
the results.


```
Rscript \
  -e '{
    library(parallel)

    QC_CUT <- 0.3
    partners <- c("HUA","REGIONH","MUW","UCY","UTARTU","HMGU","VUA","TAUH","QIMR",
        "BIB","WIS")

    mainPath <- "DIRECTORY_STRUCTURE"
    pat <- ".*[^0-9](1[0-9]|2[0-2]|[1-9])[^0-9].*\\.txt\\.gz$"
    pat2 <- ".*(1[0-9]|2[0-2]|[1-9]).*\\.txt\\.gz$"

    partnerSnps <- finalSnps <- vector("list",length(partners))
    names(partnerSnps) <- names(finalSnps) <- partners
    for (p in partners) {
        message("Reading data for partner ",p)
        if (dir.exists(file.path(mainPath,p))) {
            pops <- dir(file.path(mainPath,p),full.names=TRUE)
            pops <- pops[dir.exists(pops)]
            partnerSnps[[p]] <- lapply(pops,function(x) {
                message("  Reading dataset ",x)
                snps <- dir(x,pattern=pat,full.names=TRUE)
                if (length(snps) == 0)
                    snps <- dir(x,pattern=pat2,full.names=TRUE)
                chrSnps <- mclapply(snps,function(y) {
                    message("    Reading file ",y)
                    vars <- read.delim(y,header=FALSE,comment.char="#")
                    
                    metricsField <- strsplit(vars[,5],split=";")
                    metrics <- do.call("rbind",lapply(metricsField,function(z) {
                        tmp <- strsplit(z,split="=")
                        out <- lapply(tmp,function(a) return(as.numeric(a[2])))
                        names(out) <- lapply(tmp,function(a) return(a[1]))
                        return(unlist(out))
                    }))
                    
                    # We may have INFO, R2 or nothing
                    qcField <- ifelse("INFO" %in% colnames(metrics),"INFO",
                        ifelse("R2" %in% colnames(metrics),"R2",NA))
                    
                    if (!is.na(qcField)) {
                        keep <- metrics[,qcField] > QC_CUT
                        keepVars <- vars[keep,,drop=FALSE]
                    }
                    else
                        keepVars <- vars
                    
                    return(paste(keepVars[,1],keepVars[,2],sep=":"))
                },mc.cores=22)
                
                names(chrSnps) <- unlist(lapply(chrSnps,function(y) {
                    snp <- strsplit(y[1],":")
                    return(snp[[1]][1])
                }))
                
                return(chrSnps)
            })
            
            # We now have to apply union for each partner subset
            allChrs <- unique(unlist(lapply(partnerSnps[[p]],names)))
            # First create a list with all substudy SNPs for that chromosome
            finalSnps[[p]] <- lapply(allChrs,function(x) {
                return(unique(unlist(lapply(partnerSnps[[p]],function(y) {
                    return(y[[x]])
                }))))
            })
            names(finalSnps[[p]]) <- unlist(lapply(finalSnps[[p]],function(y) {
                snp <- strsplit(y[1],":")
                return(snp[[1]][1])
            }))
        }
    }
    nullPartner <- unlist(lapply(partnerSnps,is.null))
    finalSnps <- finalSnps[!nullPartner]

    # Now we have a list of final good SNPs per partner and per chromosome. We 
    # need to intersect good SNPs per chromosome for each partner.
    allChrs <- unique(unlist(lapply(finalSnps,names),use.names=FALSE))
    preSnps <- lapply(allChrs,function(x) {
        return(lapply(finalSnps,function(y) {
            return(y[[x]])
        }))
    })
    names(preSnps) <- allChrs

    # Final SNPs per chromosome with Reduce
    theSnps <- lapply(preSnps,function(x) {
        Reduce("intersect",x)
    })
    names(theSnps) <- allChrs

    # Write
    lapply(names(theSnps),function(n) {
        #fname <- file.path(mainPath,paste0("b4u_chr",n,"_isecqc.ids"))
        # Or if we have a standard working directory with 1000G data
        fname <- paste0("b4u_chr",n,"_isecqc.ids")
        writeLines(theSnps[[n]],fname)
    })
  }'
```

The above script yields a set of files with SNP ids named `b4u_chr*_isecqc.ids`.
We will use them to intersect with 1000 genomes SNPs so as to derive and initial
population for LD pruning.

## 5. Remove any 1000 genomes variants not present in HRC imputed cohorts

We use the list of common SNPs derived above to intersect with 1000 genomes
SNPs:

```
Rscript \
  -e '{
    # Chromosomes
    chrs <- seq_len(22)

    # Make instersections for each chromosome
    lapply(chrs,function(chr) {
        message("Finding common variants for chromosome ",chr)
        #hrcVars <- read.table(paste0("chr",chr,"_HRC_ids.txt"))[,1]
        hrcVars <- read.table(paste0("b4u_chr",chr,"_isecqc.ids"))[,1]
        kg1Vars <- read.table(paste0("ALL.chr",chr,".variant.ids"))[,1]
        commonVars <- intersect(kg1Vars,hrcVars)
        writeLines(commonVars,paste0("chr",chr,"_HRC1KG_common.ids"))
    })
  }'
```

The files `chr*_HRC1KG_common.ids` will be used with PLINK to create PLINK files
for further downstream analysis.

## 6. Convert the BCF files to PLINK format

Next, we convert BCF files to PLINK BED|BIM|FAM format so as to proceed with the
operations required for PCA:

```
for CHR in `seq 1 22`
do
  plink \
    --bcf ALL.chr"${CHR}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf \
    --keep-allele-order \
    --vcf-idspace-to _ \
    --const-fid \
    --allow-extra-chr 0 \
    --split-x b37 no-fail \
    --extract chr${CHR}_HRC1KG_common.ids \
    --make-bed \
    --out ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes
done
```

## 7. Prune variants from each chromosome

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

```
mkdir pruned

TOTAL=10000000
R2=0.95
STEP=0.05

while [ $TOTAL -gt 200000 ]
do
  R2=$(echo "$R2-$STEP" | bc)
  
    for CHR in `seq 1 22`
    do
      plink \
        --bfile ALL.chr"${CHR}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes \
        --indep-pairwise 500 50 $R2 \
        --out ./pruned/ALL.chr"${CHR}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes \
        --silent &
    done
    
    wait
  
    #TOTAL=$(wc -l ./pruned/*.in | grep total | cut -d' ' -f2)
    TOTAL=$(wc -l ./pruned/*.in | grep total | awk '{print $1}')
  
    echo "Parameters --indep-pairwise 500 50 $R2 yield $TOTAL SNPs" >> \
      pruning.log
done
```

The final value of `$TOTAL` is `0.15` and the number of SNPs is 178,135, so we
perform the final pruning:

```
for CHR in `seq 1 22`
do
  plink \
    --bfile ALL.chr"${CHR}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes \
    --indep-pairwise 500 50 0.15 \
    --out ./pruned/ALL.chr"${CHR}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes
done

for CHR in `seq 1 22`
do
  plink \
    --bfile ALL.chr"${CHR}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes \
    --extract ./pruned/ALL.chr"${CHR}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.prune.in \
    --make-bed \
    --out ./pruned/ALL.chr"${CHR}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes
done
```

## 8. Assemble the pruned SNPs and merge in preparation for PCA

Firstly, we collect a merge list file for PLINK:

```
for CHR in `seq 1 22`
do
  echo `readlink -f ./pruned/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bim` >> mergelist.txt
done
sed -i 's/\.bim//g' mergelist.txt
```

Now merge the pruned data into a single dataset for PCA using PLINK:

```
plink --merge-list mergelist.txt --out pruned_merged_1000g
```

And perform PCA using `flashpca` (the files `loads_1000g.txt` and 
`means_1000g.txt` are required for later projections of own data):

```
flashpca_x86-64 \
  --bfile pruned_merged_1000g \
  # --ndim 20 \
  --ndim 5 \
  --outpc pcs_1000g.txt \
  --outvec eig_1000g.txt \
  --outload loads_1000g.txt \
  --outmeansd means_1000g.txt
```

The files `loads_1000g.txt` and `means_1000g.txt` are required and will be used
to create projections of each partner data to the 1000 genomes PCs and these
will be used as covariates for individual partner cohort analyses. Therefore,
these files will be distributed to the partners.

Finally, create the variant list to distribute to partners

```
cut -f2 pruned_merged_1000g.bim > pca_variants.txt
```

The file `pca_variants.txt` will be distributed to the partners for the 
construction of PLINK files to be used for projection.

# Appendix I

Script provided by Anders Ericsson to extract SNP information and quality scores
from HRC imputed data:

```
#!/bin/bash
#
# Author: Anders Eriksson (anders.eriksson@ut.ee)
#
# Description: Bash script for extracting imputation and other key summary information from VCF files
# The script will generate a  (gzip compressed) files containing information from the VCF headers 
# about the meaning of the INFO field variables and, for each genetic variant in the dataset,
# the chromosome, position, reference genome allele, alternative allele and INFO field
#
# Instructions: 
# 1. Edit the script to set the PP variable to the path to the folder containing the VCF or VCF.gz with the genotype data for your cohort
# 2. Run the script: bash info_extraction_script.sh
# 3. Zip the resulting files together, label the zip file clearly with your acronym (e.g. UTARTU for Tartu) 
#    and upload to the BETTER4 Google drive folder BETTER4U > WPs > WP3 > GWAS > Code
#
# Please do not hesitate to contact Anders Eriksson (anders.eriksson@ut.ee) in case of any questions or problems
#
# Example output (from UTARTU):
# $ zcat EstBB_chr1.info.txt.gz 
# ##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated ALT Allele Frequencies">
# ##INFO=<ID=IMP,Number=0,Type=Flag,Description="Imputed marker">
# ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">
# ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
# ##INFO=<ID=INFO,Number=1,Type=Float,Description="IMPUTE2 info score">
# 1 13380   C   G   IMP;AF=0;AN=422518;AC=0;INFO=0.0220531
# 1 16071   G   A   IMP;AF=0;AN=422518;AC=0;INFO=0.130104
# 1 16141   C   T   IMP;AF=2.36676e-06;AN=422518;AC=1;INFO=0.0555347
# 1 16280   T   C   IMP;AF=2.36676e-06;AN=422518;AC=1;INFO=0.197894
# 1 49298   T   C   IMP;AF=0.766266;AN=422518;AC=323761;INFO=0.247715
# 1 54353   C   A   IMP;AF=0.000123072;AN=422518;AC=52;INFO=0.178187
# 1 54564   G   T   IMP;AF=7.10029e-06;AN=422518;AC=3;INFO=0.0864531
# 1 54591   A   G   IMP;AF=1.89341e-05;AN=422518;AC=8;INFO=0.319276
# 1 54676   C   T   IMP;AF=0.250243;AN=422518;AC=105732;INFO=0.254523

PP='path to folder containing VCF or VCF.gz files'
for FILE in `ls $PP/*.vcf $PP/*.vcf.gz`; do
   FILE2=${FILE##*/}; FILE2=${FILE2%.vcf*}
   bcftools view -h ${FILE}| grep '##INFO' | gzip -c  > ${FILE2}.info.txt.gz
   bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO\n' ${FILE} | gzip -c  >> ${FILE2}.info.txt.gz
done
```

# Open issues

No current issue.

### Notes

##### Personal script to get test variants for projection

```
cohvars <- read.delim("COHORT_imputed_merged.bim",header=FALSE)
pcavars <- read.table("./1000/pca_variants.txt")
vars <- intersect(v$V1,z$V2)
writeLines(vars,"pca_variants_sub.txt")
```

##### Commands to create DAP test 1000 genomes dataset and projection

```
plink \
  --bfile pruned_merged_1000g \
  --out for_dap_test \
  --extract ../pca_variants_sub.txt \
  --make-bed

flashpca_x86-64 \
  --bfile for_dap_test \
  --ndim 20 \
  --outpc pcs_test.txt \
  --outvec eig_test.txt \
  --outload loads_test.txt \
  --outmeansd means_test.txt

flashpca \
  --bfile COHORT_for_PCA \
  --inmeansd ./1000/means_test.txt \
  --inload ./1000/loads_test.txt \
  --project \
  --outproj proj_test.txt
```

##### Command to get a list of HRC variants from NAFLD HRC imputed data

```
for CHR in `seq 1 22`
do
  bcftools query --format '%ID' NAFLD.chr${CHR}.HRC.vcf.gz \
    > chr${CHR}_HRC_ids.txt &
done
```

##### PC projection clustering

We will need to draft an additional process for this as soon as we have some
results.
