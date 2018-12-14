# copyright 2017-2018 Regents of the University of California and the Broad Institute. All rights reserved.
#version 1.0
FROM edjuaro/pdx-hts-base:1.0

MAINTAINER Edwin Juarez <ejuarez@ucsd.edu>

ENV LANG=C LC_ALL=C

# Installing some other required libraries
USER root

RUN pip install openpyxl==2.5.12 && \
    pip install matplotlib-venn==0.11.5 && \
    pip install cuzcatlan==0.9.0 && \
    pip install ccalnoir==2.3 && \
    pip install pubchempy==1.0.4 && \
    pip install tqdm==4.28.1

# We will install these using Rscript. Not the most reproducible method, but the "R" way of doing things.
# Attempting to reproduce this:
# Installed:
# https://bioconductor.org/packages/3.6/bioc/src/contrib/affyio_1.48.0.tar.gz
# https://bioconductor.org/packages/3.6/bioc/src/contrib/annotate_1.56.2.tar.gz
# https://bioconductor.org/packages/3.6/data/annotation/src/contrib/hspeccdf_0.99.1.tar.gz
# however, I can't get the list of dependencies installed. Ugh. R...
# At least *this* container is reproducible itself
#RUN conda install -y -c bioconda bioconductor-affy==1.58.0 && \
#    conda install -y -c bioconda bioconductor-annotate

RUN mkdir /build/source
COPY install_libraries.R /build/source
RUN Rscript /build/source/install_libraries.R

USER $NB_USER

# build using this:
# docker build -t edjuaro/pdx-hts:2.0 .
