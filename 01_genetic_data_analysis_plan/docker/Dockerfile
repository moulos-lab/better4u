FROM phusion/baseimage:jammy-1.0.4

LABEL author="Stelios GKIONIS <gkionis@fleming.gr>"

# Use baseimage-docker's init system.
CMD ["/sbin/my_init"]

# Install apt-get dependencies
RUN DEBIAN_FRONTEND=noninteractive apt-get -y update && apt-get install -y -q \
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
    rm PRSice_linux.zip
ENV PATH=/tools/prsice:$PATH

WORKDIR /tools


# Install R
RUN DEBIAN_FRONTEND=noninteractive wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | \
    tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc && \
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" && \
    apt-get -y update && apt-get install -y -q r-base

# Install R dependencies
# Re organize these sections to be modular by package and to remove redundancies when not needed
# May need to add '-q', or 'dependencies=TRUE', or split command.
# May use Ncpus = 16 to speedup
RUN R -e "install.packages(c('susieR', 'parallel', 'utils', 'remotes', \
    'R.utils', 'statgenGWAS', 'rmarkdown'), Ncpus = 4)" && \
    R -e "install.packages('BiocManager')" && \
    R -e "library(BiocManager); BiocManager::install('snpStats')"

# RUN R -e "install.packages(c('ggplot2', 'data.table', 'optparse', \
#     'Rcpp', 'RcppParallel', 'RcppArmadillo', 'RcppEigen', 'Matrix', 'methods', \
#     'BH', 'SPAtest', 'SKAT', 'devtools'))"
# RUN R -e "remotes::install_github('leeshawn/MetaSKAT')"
# RUN R -e "remotes::install_github('weizhouUMICH/SAIGE')"

# Install system dependencies
RUN apt-get -y update && apt-get install -y -q \
    git pip cmake libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libcurl4-openssl-dev libssl-dev libxml2-dev libfreetype6-dev libpng-dev \
    libtiff5-dev libjpeg-dev libmariadb-dev

# # devtools is in the requirements of SAIGE
# RUN R -e "install.packages('devtools')"
# # RUN pip install cget click
# # RUN cget install statgen/savvy

# # libboost-all-dev may not be needed
# RUN apt-get -y update && apt-get install -y libpq-dev libboost-all-dev
# # # May be useless to cd
# # WORKDIR /tools/R_library
# RUN src_branch=main && \
#     repo_src_url=https://github.com/saigegit/SAIGE && \
#     git clone --depth 1 -b $src_branch $repo_src_url && \
#     Rscript ./SAIGE/extdata/install_packages.R && \
#     R CMD INSTALL SAIGE

# # Set environment variables to adjust compiler flags
# ENV CXXFLAGS="-Wno-format-security"
# # Not working
# RUN R -e "remotes::install_github('tshmak/lassosum')"

# # Copy the SAIGE executables to your container at build time
# COPY --from=wzhou88/saige:1.3.6 /usr/local/bin/step1_fitNULLGLMM.R /usr/local/bin/step1_fitNULLGLMM.R

# RUN wget https://github.com/tshmak/lassosum/releases/download/v0.4.5/lassosum_0.4.5.tar.gz
# RUN R -e "install.packages('lassosum_0.4.5.tar.gz', repos=NULL)"


# Set working directory
WORKDIR /app

# What to do by default
CMD ["bash"]

# Clean up APT when done
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


## Build instructions (add repo)
# docker build -t better4u ./
        
## Usage example
# Create the container once after download (add readonly)
# docker run -t -d --name better4u -v /path/to/project/dir:/app --network none better4u
#
# Start container if stopped
# docker start better4u
#
# Execute commands in running container
# docker exec better4u COMMAND
#
# Stop container after finishing session
# docker stop better4u

## Example commands for tools
# BCFTools
# docker exec better4u bcftools view example.vcf
#
# PRSice
# docker exec better4u PRSice.R --dir . \
#     --prsice ./PRSice \
#     --base TOY_BASE_GWAS.assoc \
#     --target TOY_TARGET_DATA \
#     --thread 1 \
#     --stat OR \
#     --binary-target T
#

# Use non root user? use mount or volume?