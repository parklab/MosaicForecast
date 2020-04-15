# largely based on rocker r-base image

FROM ubuntu:18.04

# Add user to 'staff' group, granting them write privileges to /usr/local/lib/R/site.library
RUN useradd docker \
	&& mkdir /home/docker \
	&& chown docker:docker /home/docker \
	&& addgroup docker staff

RUN apt-get update \ 
	&& apt-get install -y --no-install-recommends \
		ed \
		less \
		locales \
		vim-tiny \
		wget \
		ca-certificates \
		apt-transport-https \
		gsfonts \
		gnupg2 \
		bzip2 \
		gcc \
		git \
		libncurses-dev \
		make \
		time \
		unzip \
		vim \
		zlib1g-dev \
		liblz4-tool \
		libxt6 \
		libxml2-dev \
	&& rm -rf /var/lib/apt/lists/*

RUN apt-get -qq update && apt-get -qq -y install curl bzip2 \
    && curl -sSL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -bfp /usr/local \
    && rm -rf /tmp/miniconda.sh \
    && conda install -y python=3.6.7 \
    && conda config --add channels bioconda \
    && conda install pysam==0.15.2 -y \
    && conda install numpy==1.16.1 -y \
    && conda install pandas -y \
    && conda install pyfaidx==0.5.3 -y \
    && conda install pysamstats==1.1.2 -y \
    && conda install regex -y \
    && conda install scipy==1.2.1 -y \
    && conda update conda \
    && conda install bedtools==2.25.0 \
    && conda install vcftools \
    && apt-get -qq -y remove curl bzip2 \
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log \
    && conda clean --all --yes

RUN conda install -c rdonnellyr r-base \
    && conda install -c r r-xml \
    && conda install -c r r-ggplot2 \
    && conda install -c r r-caret \
    && conda install -c r r-randomForest \
    && conda install -c r r-e1071 \
    && conda install -c r r-glmnet \
    && conda install -c r r-rcolorbrewer \
    && conda install -c r r-devtools \
    && conda install -c r r-nnet

ENV PATH /opt/conda/bin:$PATH
RUN R -e 'install.packages("mlr", repos="http://cran.fiocruz.br/")' \ # mlr 

# set path:
ENV PATH=/usr/local/bin/samtools/:$PATH

# download other tools
WORKDIR /usr/local/bin
COPY downloads_docker.sh .
RUN . downloads_docker.sh

# 7. wrapper
COPY *.py *.R *.md ./
RUN chmod +x *py && chmod +x *.R
COPY k24.umap.wg.bw ./

 
