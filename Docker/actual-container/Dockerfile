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

USER $NB_USER

# build using this:
# docker build -t edjuaro/pdx-hts:1.0 .