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

## 2. Download 1000 Genomes PED file

```
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped ;
```

## 3. Download the GRCh37/hg19 reference genome

For better data integority, we download the human reference genome version used
by the consortium:

```
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz ;
```

## 4. Convert the 1000 Genomes files to BCF

In the following we:

* Ensure that multi-allelic calls are split and that indels are left-aligned 
compared to reference genome (1st pipe)
* Set the VCF ID field to a unique value: `CHROM:POS:REF:ALT` (2nd pipe)
* Remove duplicates (3rd pipe)

Notes:

- `-I +'%CHROM:%POS:%REF:%ALT'` means that unset IDs will be set to 
`CHROM:POS:REF:ALT*`
- `-x ID -I +'%CHROM:%POS:%REF:%ALT'` first erases the current ID and then sets 
it to `CHROM:POS:REF:ALT*`

```
for CHR in `seq 1 22`; 
do
  echo "Converting $CHR"
  
  bcftools norm --multiallelics -any --check-ref w \
    --fasta-ref human_g1k_v37.fasta \
    ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz | \
  bcftools annotate --remove ID --set-id +'%CHROM:%POS:%REF:%ALT' | \
  bcftools norm --output-type b --rm-dup both \
    --output ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf;
done
```

Then index:

```
for CHR in `seq 1 22`; 
do
  bcftools index \    
    ALL.chr"${CHR}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf;
done
```

The above command will probably take some time, so they can be put in a shell
script and executed with `nohup`.

## 5. Convert the BCF files to PLINK format

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
  --make-bed \
  --out ALL.chr"${CHR}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes
done
```

6. Prune variants from each chromosome

Prior to PCA, it is common practice to not use variants in LD, so LD pruning is
performed.

* `--maf 0.05`, only retains SNPs with MAF greater than 5%.
* `--indep [window size] [step size/variant count)] [Variance inflation factor (VIF) threshold]` controls the pruning parameters. e.g. indep `50 5 1.5`, 
generates a list of markers in approximate linkage equilibrium - takes 50 SNPs 
at a time and then shifts by 5 for the window. VIF (1/(1-r^2)) is the cut-off 
for linkage disequilibrium.

We split the commands below in two `for` loops so that `&` can be used to easily
parallelize per chromosome wherever possible:

```
mkdir pruned

for CHR in `seq 1 22`
do
  plink \
    --bfile ALL.chr"${CHR}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes \
    --maf 0.05 \
    --indep 50 5 1.5 \
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

7. Assemble the pruned SNPs and merge in preparation for PCA

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

8. Create the variant list to distribute to partners

```
cut -f2 pruned_merged_1000g.bim > pca_variants.txt
```

The file `pca_variants.txt` will be distributed to the partners for the 
construction of PLINK files to be used for projection.

# Open issues

* The above process produces >700,000 variants for PCA. They may be too many
for certain smaller cohorts, as because of size, imputation filtering may remove
a lot of variants. We could use the process described [here](https://finngen.gitbook.io/documentation/methods/phewas/quality-checks).
* This process will be executed only by the central analysis team? Probably yes,
and then a list of SNPs will be distributed.

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

##### PC projection clustering

We will need to draft an additional process for this as soon as we have some
results.
