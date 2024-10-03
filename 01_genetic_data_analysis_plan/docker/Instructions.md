Container Instructions
=============================================================================

# Download and Build Instructions

## Apptainer/Singularity
Download and build an Apptainer image from the Docker image on Docker Hub.
For convenience put it in the project dir.
```bash
apptainer pull docker://stgkionis/better4u 
```

## Docker
Download from Docker Hub.
```bash
docker pull docker://stgkionis/better4u 
```

# Usage Instructions

## With Apptainer/Singularity
### Create the Container
Create the container with no external network connectivity with `--net --network none`.
Do this from the project dir so it is mounted automatically.
```bash
apptainer run --net --network none ./better4u-latest.sif
```

### Execute Commands in the Running Container
Simply run commands as you would in any interactive shell and `exit` when finished.

## With Docker
### Create the Container
Create the container with no external network connectivity with `--network none`
and mount the project dir.
```bash
docker run -t -d --name better4u -v /path/to/project/dir:/app --network none better4u
```

### Start the Container
If the container is stopped, start it again.
```bash
docker start better4u
```

### Execute Commands in the Running Container
See tool-specific examples below:
```bash
docker exec better4u <COMMAND>
```

### Stop the Container
Stop the container when you're done.
```bash
docker stop better4u
```

# Example Commands for Tools

Commands are given for Apptainer/Singularity. For use with Docker, simply prepend `docker exec better4u <COMMAND>`.
### Example commands
```bash
bcftools --help
vcftools --help
# Do not change `--prsice`, add additional arguments like `--help` afterwards:
PRSice.R --prsice /tools/prsice/PRSice_linux --help
plink --help
regenie --help
gcta64
snptest -help
flashpca --help
MR-MEGA --help
# Use the utility scripts like this (per documentation):
R --slave --vanilla < /tools/mrmega/fixP.r
```

## Tool Versions
Here are the versions of the tools included in this Docker image:
| Tool       | Version |
|------------|---------|
| bcftools   | 1.20    |
| htslib     | 1.20    |
| VCFtools   | 0.1.16  |
| PLINK      | 1.9     |
| REGENIE    | 3.5     |
| GCTA       | 1.94.1  |
| SNPTEST    | 2.5.6   |
| PRSice2    | 2.3.5   |
| MR-MEGA    | 0.2     |
| FlashPCA   | 2.0     |
| R          | 4.4.1   |

## R Packages
### User-Specified Packages
Here are the versions of the user-specified R packages:
| Package      | Version |
|--------------|---------|
| susieR       | 0.12.35 |
| snpStats     | 1.54.0  |
| parallel     | 4.4.1   |
| utils        | 4.4.1   |
| remotes      | 2.5.0   |
| R.utils      | 2.12.3  |
| statgenGWAS  | 1.0.9   |
| rmarkdown    | 2.27    |
| SAIGE        | 1.3.6   |
| lassosum     | 0.4.5   |
