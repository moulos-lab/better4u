Calculation of main Principal Components for BETTER4U
================================================================================

The following is based on [this](https://www.biostars.org/p/335605/) tutorial
from [Kevin Blighe](https://github.com/kevinblighe) in Biostars.

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

## 4. Remove any 1000 genomes variants not present in HRC

Firstly, retrieve a list of all the HRC panel variants from [here](#). Then for
each chromosome we produce a list of common variants (should be most of them):

```
Rscript \
  -e '{
    # Chromosomes
    chrs <- seq_len(22)

    # Make instersections for each chromosome
    lapply(chrs,function(chr) {
        message("Finding common variants for chromosome ",chr)
        hrcVars <- read.table(paste0("chr",chr,"_HRC_ids.txt"))[,1]
        kg1Vars <- read.table(paste0("ALL.chr",chr,".variant.ids"))[,1]
        commonVars <- intersect(kg1Vars,hrcVars)
        writeLines(commonVars,paste0("chr",chr,"_HRC1KG_common.ids"))
    })
  }'
```

The files `chr*_HRC1KG_common.ids` will be used with PLINK to create PLINK files
for further downstream analysis.

## 4. Convert the BCF files to PLINK format

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

## 5. Prune variants from each chromosome

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

## 6. Assemble the pruned SNPs and merge in preparation for PCA

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
  --ndim 20 \
  --outpc pcs_1000g.txt \
  --outvec eig_1000g.txt \
  --outload loads_1000g.txt \
  --outmeansd means_1000g.txt
```

The files `loads_1000g.txt` and `means_1000g.txt` are required and will be used
to create projections of each partner data to the 1000 genomes PCs and these
will be used as covariates for individual partner cohort analyses. Therefore,
these files will be distributed to the partners.

7. Create the variant list to distribute to partners

```
cut -f2 pruned_merged_1000g.bim > pca_variants.txt
```

The file `pca_variants.txt` will be distributed to the partners for the 
construction of PLINK files to be used for projection.

# Open issues

* We need to make sure that the set of variants for PCA projection exists in all
cohorts. To this end, step 4 could be expanded to use the intersections of high
quality HRC variants from each partner which we need to receive (so far only
from UTARTU - we haven't asked yet).

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
