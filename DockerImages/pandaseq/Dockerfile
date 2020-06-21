FROM ubuntu:14.04

MAINTAINER Per Greisen <pjug@novonordisk.com>

# Initial configuration
RUN apt-get update --fix-missing && apt-get install -y build-essential wget \
    bzip2 ca-certificates libglib2.0-0 libxext6 libsm6 libxrender1 pigz s3cmd \
    git mercurial subversion libtool automake zlib1g-dev libbz2-dev pkg-config \
    muscle mafft cd-hit unzip libfontconfig1
RUN ln -s /usr/bin/cdhit /usr/bin/cd-hit
RUN mkdir /tools

# Anaconda
RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/archive/Anaconda2-4.0.0-Linux-x86_64.sh && \
    /bin/bash /Anaconda2-4.0.0-Linux-x86_64.sh -b -p /opt/conda && \
    rm /Anaconda2-4.0.0-Linux-x86_64.sh
ENV PATH /opt/conda/bin:$PATH

# PANDAseq
RUN cd /tools && \
    git clone https://github.com/neufeld/pandaseq && \
    cd pandaseq && \
    ./autogen.sh && \
    ./configure && \
    make && \
    make install && \
    ldconfig

# BaseSpace Python SDK
RUN cd /tools && \
    git clone https://github.com/basespace/basespace-python-sdk && \
    cd basespace-python-sdk/src && \
    python setup.py install

# FASTQC
RUN cd /tools && \
    wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip && \
    unzip fastqc_v0.11.5.zip && \
    ln -s FastQC/fastqc /usr/local/bin/fastqc

# Cutadapt
RUN pip install --upgrade pip
RUN pip install cutadapt
