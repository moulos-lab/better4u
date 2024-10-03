FROM phusion/baseimage:jammy-1.0.4 AS base

LABEL description="Container with tools for the Better4u project"
LABEL maintainer="Stelios GKIONIS <gkionis@fleming.gr>"

# # baseimage-docker's init system cannot be used with Apptainer.
# ENTRYPOINT ["/sbin/my_init"]

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

# PRSice installation differs from github instructions
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

# Install SAIGE and dependencies
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

# Install MR-MEGA
WORKDIR /tools/mrmega

RUN wget https://tools.gi.ut.ee/tools/MR-MEGA_v0.2.zip && \
    unzip MR-MEGA_v0.2.zip && \
    rm MR-MEGA_v0.2.zip && \
    make && \
    curl -o fixP.r https://tools.gi.ut.ee/tools/fixP.r && \
    curl -o manh.r https://tools.gi.ut.ee/tools/manh.r && \
    curl -o qq.r https://tools.gi.ut.ee/tools/qq.r && \
    sed -i '1i#!/usr/bin/Rscript --slave --vanilla\n' fixP.r manh.r qq.r && \
    chmod +x fixP.r manh.r qq.r
ENV PATH=/tools/mrmega:$PATH

# Install FlashPCA
WORKDIR /tools

RUN apt-get update && apt-get install -y -q \
    libeigen3-dev libboost-program-options-dev && \
    wget https://github.com/yixuan/spectra/archive/refs/tags/v0.9.0.tar.gz && \
    tar -xzvf v0.9.0.tar.gz && rm v0.9.0.tar.gz && \
    wget https://github.com/gabraham/flashpca/releases/download/v2.0/flashpca_2.0.tar.gz && \
    tar -xzvf flashpca_2.0.tar.gz && rm flashpca_2.0.tar.gz && \
    cd flashpca_2.0 && \
    sed -i '/^CXXFLAGS =/s/$/ -fpermissive/' Makefile && \
    make all EIGEN_INC=/usr/include/eigen3 \
        BOOST_INC=/usr/include/boost \
        BOOST_LIB=/usr/lib/x86_64-linux-gnu \
        SPECTRA_INC=/tools/spectra-0.9.0/include/Spectra
ENV PATH=/tools/flashpca_2.0:$PATH

# Install lassosum
RUN sed -i 's/-Werror=format-security/-Wno-format-security/g' /etc/R/Makeconf && \
    R -e "library(remotes); install_github('tshmak/lassosum@v0.4.5')" && \
    sed -i 's/-Wno-format-security/-Werror=format-security/g' /etc/R/Makeconf

# Reset frontend to default
ENV DEBIAN_FRONTEND=teletype

# Set working directory for Docker
WORKDIR /app

# What to do by default
CMD ["bash"]

# Clean up APT when done
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
