FROM phusion/baseimage:jammy-1.0.4 AS base

LABEL author="Stelios GKIONIS <gkionis@fleming.gr>"

# Use baseimage-docker's init system.
CMD ["/sbin/my_init"]

# Make frontend noninteractive for setup
ENV DEBIAN_FRONTEND=noninteractive

# Install apt-get dependencies for command line tools
RUN apt-get update && apt-get install -y -q \
    gcc zlib1g-dev libbz2-dev liblzma-dev g++ make pkg-config \
    wget zip unzip lbzip2 software-properties-common dirmngr

# Install and configure command line tools
WORKDIR /tools

RUN wget https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2 && \
    tar -xf bcftools-1.20.tar.bz2 && \
    rm bcftools-1.20.tar.bz2 && \
    cd bcftools-1.20 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && \
    rm -r bcftools-1.20

RUN wget https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2 && \
    tar -xf htslib-1.20.tar.bz2 && \
    rm htslib-1.20.tar.bz2 && \
    cd ./htslib-1.20 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && \
    rm -r htslib-1.20

RUN PERL5LIB=/tools/vcftools-0.1.16/src/perl/ \
    wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz && \
    tar -xf vcftools-0.1.16.tar.gz && \
    rm vcftools-0.1.16.tar.gz && \
    cd ./vcftools-0.1.16 && \
    ./configure && \
    make && \
    make install

WORKDIR /tools/plink

RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip && \
    unzip plink_linux_x86_64_20231211.zip && \
    rm plink_linux_x86_64_20231211.zip 
ENV PATH=/tools/plink:$PATH

WORKDIR /tools/regenie

RUN wget https://github.com/rgcgithub/regenie/releases/download/v3.5/regenie_v3.5.gz_x86_64_Linux.zip && \
    unzip regenie_v3.5.gz_x86_64_Linux.zip && \
    rm regenie_v3.5.gz_x86_64_Linux.zip && \
    mv regenie_v3.5.gz_x86_64_Linux regenie
ENV PATH=/tools/regenie:$PATH

WORKDIR /tools

RUN wget https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.1-linux-kernel-3-x86_64.zip && \
    unzip gcta-1.94.1-linux-kernel-3-x86_64.zip && \
    rm gcta-1.94.1-linux-kernel-3-x86_64.zip && \
    rm -r __MACOSX && \
    mv gcta-1.94.1-linux-kernel-3-x86_64 gcta
ENV PATH=/tools/gcta:$PATH

RUN wget http://www.well.ox.ac.uk/~gav/resources/snptest_v2.5.6_CentOS_Linux7.8-x86_64_dynamic.tgz && \
    tar -xf snptest_v2.5.6_CentOS_Linux7.8-x86_64_dynamic.tgz && \
    rm snptest_v2.5.6_CentOS_Linux7.8-x86_64_dynamic.tgz && \
    mv snptest_v2.5.6_CentOS_Linux7.8.2003-x86_64_dynamic snptest && \
    mv snptest/snptest_v2.5.6 snptest/snptest
ENV PATH=/tools/snptest:$PATH

# Differs from github instructions
WORKDIR /tools/prsice

RUN wget https://github.com/choishingwan/PRSice/releases/download/2.3.5/PRSice_linux.zip && \
    unzip PRSice_linux.zip && \
    rm PRSice_linux.zip && \
    chmod +x /tools/prsice/PRSice_linux
ENV PATH=/tools/prsice:$PATH

WORKDIR /tools

# Install R
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | \
    tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc && \
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" && \
    apt-get -y update && apt-get install -y -q r-base

# Install R packages 
RUN R -e "install.packages(c('susieR', 'parallel', 'utils', 'remotes', \
    'R.utils', 'statgenGWAS', 'rmarkdown'), Ncpus = 4)" && \
    R -e "install.packages('BiocManager')" && \
    R -e "library(BiocManager); BiocManager::install('snpStats')"

# Install dependencies and SAIGE
RUN apt-get update && apt-get install -y -q \
    cmake git pip file libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libcurl4-openssl-dev libssl-dev libxml2-dev libfreetype6-dev libpng-dev \
    libtiff5-dev libjpeg-dev libmariadb-dev libpq-dev && \
    src_branch=main && \
    repo_src_url=https://github.com/saigegit/SAIGE && \
    git clone --depth 1 -b $src_branch $repo_src_url && \
    Rscript ./SAIGE/extdata/install_packages.R && \
    # Manually run pip because SAIGE configure file sets -eu and fails otherwise \
    # Alternatively set -- to unset -eu as it only fails at pip warning \
    pip install cget click && \
    R CMD INSTALL SAIGE && \
    # RUN R CMD INSTALL SAIGE --configure-args='SHELL_OPTS="set --"' \
    chmod +x /tools/SAIGE/extdata/step* && \
    cp /tools/SAIGE/extdata/step* /usr/local/bin/

# Reset frontend to default
ENV DEBIAN_FRONTEND=teletype

# Set working directory
WORKDIR /app

# # What to do by default
# CMD ["bash"]

# Clean up APT when done
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


# =============================================================================
# Docker Container Instructions
# =============================================================================

# ## Build Instructions
# Build the Docker image with a custom tag
# ```bash
# docker build -t better4u ./
# ```

# ## Usage Example
# ### Create the Container
# Run the container with no external network connectivity (`--network none`)
# ```bash
# docker run -t -d --name better4u -v /path/to/project/dir:/app --network none better4u
# ```

# ### Start the Container
# If the container is stopped, start it again
# ```bash
# docker start better4u
# ```

# ### Execute Commands in the Running Container
# To execute commands inside the running container (see examples below):
# ```bash
# docker exec better4u <COMMAND>
# ```

# ### Stop the Container
# Stop the container when you're done
# ```bash
# docker stop better4u
# ```

# ## Example Commands for Tools
# ### BCFTools
# ```bash
# docker exec better4u bcftools --help
# ```

# ### VCFtools
# ```bash
# docker exec better4u vcftools --help
# ```

# ### PRSice2
# Do not change the first line add additional arguments like `--help`
# ```bash
# docker exec better4u PRSice.R --prsice /tools/prsice/PRSice_linux \
#     --help
# ```

# ### PLINK
# ```bash
# docker exec better4u plink --help
# ```

# ### REGENIE
# ```bash
# docker exec better4u regenie --help
# ```

# ### GCTA
# ```bash
# docker exec better4u gcta64
# ```

# ### SNPTEST
# ```bash
# docker exec better4u snptest -help
# ```

# ## Tool Versions
# Here are the versions of the tools included in this Docker image:
# ```
# | Tool       | Version |
# |------------|---------|
# | bcftools   | 1.20    |
# | htslib     | 1.20    |
# | VCFtools   | 0.1.16  |
# | PLINK      | 1.9     |
# | REGENIE    | 3.5     |
# | GCTA       | 1.94.1  |
# | SNPTEST    | 2.5.6   |
# | PRSice2    | 2.3.5   |
# | R          | 4.4.1   |
# ```

# ## R Packages
# ### User-Specified Packages
# Here are the versions of the user-specified R packages:
# ```
# | Package      | Version |
# |--------------|---------|
# | susieR       | 0.12.35 |
# | snpStats     | 1.54.0  |
# | parallel     | 4.4.1   |
# | utils        | 4.4.1   |
# | remotes      | 2.5.0   |
# | R.utils      | 2.12.3  |
# | statgenGWAS  | 1.0.9   |
# | rmarkdown    | 2.27    |
# | SAIGE        | 1.3.6   |
# ```

# ### All Installed Packages
# Here are the versions of all R packages installed including dependencies:
# ```
# | Package            | Version  |
# |--------------------|----------|
# | askpass            | 1.2.0    |
# | backports          | 1.5.0    |
# | base64enc          | 0.1-3    |
# | bench              | 1.1.3    |
# | BH                 | 1.84.0-0 |
# | BiocGenerics       | 0.50.0   |
# | BiocManager        | 1.30.23  |
# | BiocVersion        | 3.19.1   |
# | bit                | 4.0.5    |
# | bit64              | 4.0.5    |
# | blob               | 1.2.4    |
# | brew               | 1.0-10   |
# | brio               | 1.1.5    |
# | broom              | 1.0.6    |
# | bslib              | 0.7.0    |
# | cachem             | 1.1.0    |
# | callr              | 3.7.6    |
# | cli                | 3.6.3    |
# | clipr              | 0.8.0    |
# | colorspace         | 2.1-0    |
# | commonmark         | 1.9.1    |
# | covr               | 3.6.4    |
# | cpp11              | 0.4.7    |
# | crayon             | 1.5.3    |
# | credentials        | 2.0.1    |
# | crosstalk          | 1.2.1    |
# | curl               | 5.2.1    |
# | cyclocomp          | 1.1.1    |
# | data.table         | 1.15.4   |
# | DBI                | 1.2.3    |
# | dbplyr             | 2.5.0    |
# | desc               | 1.4.3    |
# | devtools           | 2.4.5    |
# | diffobj            | 0.3.5    |
# | digest             | 0.6.36   |
# | docopt             | 0.7.1    |
# | downlit            | 0.4.4    |
# | dplyr              | 1.1.4    |
# | DT                 | 0.33     |
# | ellipsis           | 0.3.2    |
# | evaluate           | 0.24.0   |
# | fansi              | 1.0.6    |
# | farver             | 2.1.2    |
# | fastmap            | 1.2.0    |
# | foghorn            | 1.6.0    |
# | fontawesome        | 0.5.2    |
# | fs                 | 1.6.4    |
# | gargle             | 1.5.2    |
# | generics           | 0.1.3    |
# | gert               | 2.1.0    |
# | getopt             | 1.20.4   |
# | ggplot2            | 3.5.1    |
# | gh                 | 1.4.1    |
# | gitcreds           | 0.1.2    |
# | glue               | 1.7.0    |
# | gmailr             | 2.0.0    |
# | gtable             | 0.3.5    |
# | here               | 1.0.1    |
# | highr              | 0.11     |
# | htmltools          | 0.5.8.1  |
# | htmlwidgets        | 1.6.4    |
# | httpuv             | 1.6.15   |
# | httr               | 1.4.7    |
# | httr2              | 1.0.2    |
# | hunspell           | 3.0.3    |
# | ini                | 0.3.1    |
# | inline             | 0.3.19   |
# | irlba              | 2.3.5.1  |
# | isoband            | 0.2.7    |
# | jquerylib          | 0.1.4    |
# | jsonlite           | 1.8.8    |
# | knitr              | 1.48     |
# | labeling           | 0.4.3    |
# | Lahman             | 11.0-0   |
# | later              | 1.3.2    |
# | lazyeval           | 0.2.2    |
# | lifecycle          | 1.0.4    |
# | lintr              | 3.1.2    |
# | lobstr             | 1.1.2    |
# | magrittr           | 2.0.3    |
# | markdown           | 1.13     |
# | Matrix             | 1.7-0    |
# | matrixStats        | 1.3.0    |
# | memoise            | 2.0.1    |
# | MetaSKAT           | 0.81     |
# | microbenchmark     | 1.4.10   |
# | mime               | 0.12     |
# | miniUI             | 0.1.1.1  |
# | mixsqp             | 0.3-54   |
# | mockery            | 0.4.4    |
# | munsell            | 0.5.1    |
# | nycflights13       | 1.0.2    |
# | openssl            | 2.2.0    |
# | optparse           | 1.7.5    |
# | pillar             | 1.9.0    |
# | pingr              | 2.0.3    |
# | pkgbuild           | 1.4.4    |
# | pkgconfig          | 2.0.3    |
# | pkgdown            | 2.1.0    |
# | pkgKitten          | 0.2.3    |
# | pkgload            | 1.4.0    |
# | plogr              | 0.2.0    |
# | plyr               | 1.8.9    |
# | png                | 0.1-8    |
# | praise             | 1.0.0    |
# | prettyunits        | 1.2.0    |
# | processx           | 3.8.4    |
# | profmem            | 0.6.0    |
# | profvis            | 0.3.8    |
# | promises           | 1.3.0    |
# | ps                 | 1.7.7    |
# | purrr              | 1.0.2    |
# | qlcMatrix          | 0.9.8    |
# | R.methodsS3        | 1.8.2    |
# | R.oo               | 1.26.0   |
# | R.utils            | 2.12.3   |
# | R6                 | 2.5.1    |
# | ragg               | 1.3.2    |
# | rappdirs           | 0.3.3    |
# | rcmdcheck          | 1.4.0    |
# | RColorBrewer       | 1.1-3    |
# | Rcpp               | 1.0.13   |
# | RcppArmadillo      | 14.0.0-1 |
# | RcppEigen          | 0.3.4.0.0|
# | RcppParallel       | 5.1.8    |
# | RcppProgress       | 0.4.2    |
# | RcppTOML           | 0.2.2    |
# | rematch            | 2.0.0    |
# | rematch2           | 2.1.2    |
# | remotes            | 2.5.0    |
# | reshape            | 0.8.9    |
# | reticulate         | 1.38.0   |
# | rex                | 1.2.1    |
# | RhpcBLASctl        | 0.23-42  |
# | rhub               | 2.0.0    |
# | rlang              | 1.1.4    |
# | rmarkdown          | 2.27     |
# | RMySQL             | 0.10.27  |
# | roxygen2           | 7.3.2    |
# | RPostgreSQL        | 0.7-6    |
# | rprojroot          | 2.0.4    |
# | RSpectra           | 0.16-2   |
# | RSQLite            | 2.3.7    |
# | rstudioapi         | 0.16.0   |
# | RUnit              | 0.4.33   |
# | rversions          | 2.1.2    |
# | rvest              | 1.0.4    |
# | SAIGE              | 1.3.6    |
# | sass               | 0.4.9    |
# | scales             | 1.3.0    |
# | selectr            | 0.4-2    |
# | sessioninfo        | 1.2.2    |
# | sfsmisc            | 1.1-18   |
# | shiny              | 1.8.1.1  |
# | SKAT               | 2.2.5    |
# | slam               | 0.1-51   |
# | snpStats           | 1.54.0   |
# | sommer             | 4.3.4    |
# | sourcetools        | 0.1.7-1  |
# | sparsesvd          | 0.2-2    |
# | SPAtest            | 3.1.2    |
# | spelling           | 2.3.0    |
# | statgenGWAS        | 1.0.9    |
# | stringi            | 1.8.4    |
# | stringr            | 1.5.1    |
# | susieR             | 0.12.35  |
# | sys                | 3.4.2    |
# | systemfonts        | 1.1.0    |
# | testthat           | 3.2.1.1  |
# | textshaping        | 0.4.0    |
# | tibble             | 3.2.1    |
# | tidyr              | 1.3.1    |
# | tidyselect         | 1.2.1    |
# | tinytest           | 1.4.1    |
# | tinytex            | 0.52     |
# | urlchecker         | 1.0.1    |
# | usethis            | 2.2.3    |
# | utf8               | 1.2.4    |
# | vctrs              | 0.6.5    |
# ```

# =============================================================================