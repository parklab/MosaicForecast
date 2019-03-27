FROM ubuntu:16.04
# 1. general updates & installing necessary Linux components
RUN apt-get update -y && apt-get install -y \
    bzip2 \
    gcc \
    git \
    less \
    libncurses-dev \
    make \
    time \
    unzip \
    vim \
    wget \
    zlib1g-dev \
    liblz4-tool

# 2. conda and pysam
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda3 -b
ENV PATH=/miniconda3/bin:$PATH
RUN conda update -y conda \
    && rm Miniconda3-latest-Linux-x86_64.sh
RUN conda config --add channels r \
    && conda config --add channels bioconda \
    && conda install pysam==0.11.2.2 -y \
    && conda install numpy==1.16.1 -y \
    && conda install pyfaidx==0.5.3 -y \
    && conda install pysamstats==1.0.1 -y \
    && conda install scipy==1.2.1 -y 

## 3. R packages:
## http://o2r.info/containerit/articles/containerit.html

FROM rocker/r-ver:3.4.1
RUN export DEBIAN_FRONTEND=noninteractive; apt-get -y update \
  && apt-get install -y libxml2-dev
RUN ["install2.r", "assertthat", "backports", "BBmisc", "bindr", "bindrcpp", "broom", "caret", "checkmate", "class", "codetools", "colorspace", "CVST", "data.table", "ddalpha", "DEoptimR", "dimRed", "dplyr", "DRR", "e1071", "fastmatch", "foreach", "foreign", "ggplot2", "glmnet", "glue", "gower", "gtable", "ipred", "iterators", "kernlab", "lattice", "lava", "lazyeval", "lubridate", "magrittr", "MASS", "Matrix", "mlr", "mnormt", "ModelMetrics", "munsell", "nlme", "nnet", "parallelMap", "ParamHelpers", "pillar", "pkgconfig", "plyr", "prodlim", "psych", "purrr", "R6", "RColorBrewer", "Rcpp", "RcppRoll", "recipes", "reshape2", "rlang", "robustbase", "rpart", "scales", "sfsmisc", "stringi", "stringr", "survival", "tibble", "tidyr", "tidyselect", "timeDate", "withr", "XML", "randomForest"]
#WORKDIR /payload/
#CMD ["R", "--vanilla", "-f", "test.R"]

#FROM rocker/r-ver:3.4.1
#RUN ["install2.r", "assertthat", "bindr", "bindrcpp", "broom", "caret", "class", "codetools", "colorspace", "CVST", "ddalpha", "DEoptimR", "dimRed", "dplyr", "DRR", "e1071", "foreach", "foreign", "ggplot2", "glue", "gower", "gtable", "ipred", "iterators", "kernlab", "lattice", "lava", "lazyeval", "lubridate", "magrittr", "MASS", "Matrix", "mnormt", "ModelMetrics", "munsell", "nlme", "nnet", "pillar", "pkgconfig", "plyr", "prodlim", "psych", "purrr", "R6", "randomForest", "Rcpp", "RcppRoll", "recipes", "reshape2", "rlang", "robustbase", "rpart", "scales", "sfsmisc", "stringi", "stringr", "survival", "tibble", "tidyr", "tidyselect", "timeDate", "withr"]
#WORKDIR /payload/
#CMD ["R", "--vanilla", "-f", "test.R"]


# 4. download other tools
WORKDIR /usr/local/bin
COPY downloads.sh .
RUN . downloads.sh

# 5. set path
ENV PATH=/usr/local/bin/bwa/:$PATH
ENV PATH=/usr/local/bin/samtools/:$PATH

# 6. supporting UTF-8
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

# 7. wrapper
COPY *.py *.R *.md ./
RUN chmod +x *py
#COPY *.pyc *.sh ./
#RUN chmod +x *.pyc

# 8. default command
CMD ["ls /usr/local/bin"]

