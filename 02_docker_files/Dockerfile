FROM phusion/baseimage:jammy-1.0.4

LABEL mainterer.name="Stelios GKIONIS"
LABEL maintener.email="gkionis@fleming.gr"

# Use baseimage-docker's init system.
CMD ["/sbin/my_init"]

RUN apt-get -y update

# Install apt-get dependencies
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y -q \
    gcc zlib1g-dev libbz2-dev liblzma-dev g++ make pkg-config \
    wget zip unzip lbzip2 software-properties-common dirmngr


# Install and configure command line tools
WORKDIR /tools

RUN wget https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2
RUN tar -xf bcftools-1.20.tar.bz2
WORKDIR ./bcftools-1.20
RUN ./configure --prefix=/usr/local
RUN make
RUN make install

WORKDIR /tools
RUN rm -r bcftools-1.20

RUN wget https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2
RUN tar -xf htslib-1.20.tar.bz2
RUN rm htslib-1.20.tar.bz2
WORKDIR ./htslib-1.20
RUN ./configure --prefix=/usr/local
RUN make
RUN make install

WORKDIR /tools
RUN rm -r htslib-1.20

RUN wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz
RUN tar -xf vcftools-0.1.16.tar.gz
RUN rm vcftools-0.1.16.tar.gz
ENV PERL5LIB=/tools/vcftools-0.1.16/src/perl/
WORKDIR ./vcftools-0.1.16
RUN ./configure
RUN make
RUN make install

WORKDIR /tools/plink

RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip
RUN unzip plink_linux_x86_64_20231211.zip
RUN rm plink_linux_x86_64_20231211.zip
ENV PATH=/tools/plink:$PATH

WORKDIR /tools/regenie

RUN wget https://github.com/rgcgithub/regenie/releases/download/v3.5/regenie_v3.5.gz_x86_64_Linux.zip
RUN unzip regenie_v3.5.gz_x86_64_Linux.zip
RUN rm regenie_v3.5.gz_x86_64_Linux.zip
RUN mv regenie_v3.5.gz_x86_64_Linux regenie
ENV PATH=/tools/regenie:$PATH

WORKDIR /tools

RUN wget https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.1-linux-kernel-3-x86_64.zip
RUN unzip gcta-1.94.1-linux-kernel-3-x86_64.zip
RUN rm gcta-1.94.1-linux-kernel-3-x86_64.zip
RUN rm -r __MACOSX
RUN mv gcta-1.94.1-linux-kernel-3-x86_64 gcta
ENV PATH=/tools/gcta:$PATH

RUN wget http://www.well.ox.ac.uk/~gav/resources/snptest_v2.5.6_CentOS_Linux7.8-x86_64_dynamic.tgz
RUN tar -xf snptest_v2.5.6_CentOS_Linux7.8-x86_64_dynamic.tgz
RUN rm snptest_v2.5.6_CentOS_Linux7.8-x86_64_dynamic.tgz
RUN mv snptest_v2.5.6_CentOS_Linux7.8.2003-x86_64_dynamic snptest
RUN mv snptest/snptest_v2.5.6 snptest/snptest
ENV PATH=/tools/snptest:$PATH

WORKDIR /tools/prsice

RUN wget https://github.com/choishingwan/PRSice/releases/download/2.3.5/PRSice_linux.zip
RUN unzip PRSice_linux.zip
RUN rm PRSice_linux.zip
ENV PATH=/tools/prsice:$PATH

WORKDIR /tools


# Install R
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | \
    tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc

RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"

RUN DEBIAN_FRONTEND=noninteractive apt-get install -y -q r-base

# Install R dependencies
# May need to add '-q', or 'dependencies=TRUE', or split command.
RUN R -e "install.packages(c('susieR', 'parallel', 'utils', 'remotes', \
    'R.utils', 'statgenGWAS', 'rmarkdown', 'ggplot2', 'data.table', 'optparse', \
    'Rcpp', 'RcppParallel', 'RcppArmadillo', 'RcppEigen', 'Matrix', 'methods', \
    'BH', 'SPAtest', 'SKAT', 'devtools'))"
RUN R -e "install.packages('BiocManager')" && \
    R -e "library(BiocManager); BiocManager::install('snpStats')"

RUN R -e "remotes::install_github('leeshawn/MetaSKAT')"
RUN R -e "remotes::install_github('weizhouUMICH/SAIGE')"
RUN R -e "remotes::install_github('tshmak/lassosum')"

# # Copy the SAIGE executables to your container at build time
# COPY --from=wzhou88/saige:1.3.6 /usr/local/bin/step1_fitNULLGLMM.R /usr/local/bin/step1_fitNULLGLMM.R

# RUN wget https://github.com/tshmak/lassosum/releases/download/v0.4.5/lassosum_0.4.5.tar.gz
# RUN R -e "install.packages('lassosum_0.4.5.tar.gz', repos=NULL)"


# Set working directory
WORKDIR /app

# What to do by default
CMD ["bash"]



# PRSice.R --prsice /tools/prsice/PRSice_linux

# Clean up APT when done
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*