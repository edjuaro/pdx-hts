# copyright 2017-2018 Regents of the University of California and the Broad Institute. All rights reserved.
#version 2.0
FROM jupyter/datascience-notebook:1085ca054a5f

MAINTAINER Edwin Juarez <ejuarez@ucsd.edu>

ENV LANG=C LC_ALL=C

# Installing some required libraries
RUN conda install -y rpy2==2.9.3 && \
    conda install -y imageio==2.3.0 && \
    conda install -y python-pptx==0.6.8 && \
    conda install -y -c bioconda bioconductor-biomart==2.34.2 && \
    conda install -y -c bioconda bioconductor-rhdf5==2.22.0 && \
    conda install -y -c r r-devtools==1.13.2 && \
    conda install -y -c bioconda r-sleuth==0.29.0 && \
    conda install -y -c conda-forge pytables==3.4.3 && \
    conda install -y --no-deps -c anaconda joblib==0.11 && \
    pip install pyqt5==5.10.1

# Adding GP Notebook extension

USER root
# Install the GenePattern Notebook environment and Enable extension
RUN pip install nbtools genepattern-python genepattern-notebook jupyter-wysiwyg && \
    jupyter nbextension enable --py widgetsnbextension && \
    jupyter nbextension install --py nbtools && \
    jupyter nbextension enable --py nbtools && \
    jupyter serverextension enable --py nbtools && \
    jupyter nbextension install --py genepattern && \
    jupyter nbextension enable --py genepattern && \
    jupyter serverextension enable --py genepattern && \
    jupyter nbextension install --py jupyter_wysiwyg && \
    jupyter nbextension enable --py jupyter_wysiwyg

# Installing collapsible_headings (et al.)
RUN conda install -y -c conda-forge jupyter_contrib_nbextensions && \
    jupyter nbextension enable collapsible_headings/main --user

# 2018-09-05 adding Kallisto et al.
ADD to_add /build/

# Enabling execute time
RUN jupyter nbextension enable execute_time/ExecuteTime --user

USER $NB_USER

# build using this:
# docker build -t edjuaro/pdx-hts-base:2.0 .
