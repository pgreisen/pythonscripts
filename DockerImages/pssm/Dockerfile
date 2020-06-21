FROM ubuntu:latest
# Contact me for help!
MAINTAINER pjug@novonordisk.com

ARG DEBIAN_FRONTEND=noninteractive

RUN echo BIOCONDUCTOR_DOCKER_VERSION=$BIOCONDUCTOR_DOCKER_VERSION >> /etc/environment \
	&& echo BIOCONDUCTOR_DOCKER_VERSION=$BIOCONDUCTOR_DOCKER_VERSION >> /root/.bashrc

# nuke cache dirs before installing pkgs; tip from Dirk E fixes broken img
RUN rm -f /var/lib/dpkg/available && rm -rf  /var/cache/apt/*

# issues with '/var/lib/dpkg/available' not found
# this will recreate
RUN dpkg --clear-avail

# This is to avoid the error
# 'debconf: unable to initialize frontend: Dialog'
ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends apt-utils
RUN apt-get install -y build-essential;
RUN apt-get install -y emacs wget;
RUN apt-get install -y zlib1g-dev libgl1-mesa-dev;
RUN apt-get install -y awscli;
RUN apt-get install -y zip ;
RUN apt-get install -y r-base-core;


# Install anaconda
RUN wget https://repo.continuum.io/archive/Anaconda3-2020.02-Linux-x86_64.sh

RUN bash Anaconda3-2020.02-Linux-x86_64.sh -b -p /home/anaconda
RUN rm Anaconda3-2020.02-Linux-x86_64.sh
ENV PATH /home/anaconda/bin:$PATH

COPY requirements.txt /tmp/requirements.txt
RUN pip install --requirement /tmp/requirements.txt

WORKDIR /home/

COPY data.table /home/
RUN Rscript -e "install.packages('data.table')"
# inserted to build directly into the docker images 
RUN conda install -c bioconda bioconductor-biocinstaller
RUN conda install -c conda-forge r-biocmanager

RUN  R -e "BiocManager::install('Logolas', update=FALSE, ask=FALSE)"


COPY Visualize_PSSM_w_matrix.r /home/Visualize_PSSM_w_matrix.r

COPY PSSM_analysis.py /home/PSSM_analysis.py

COPY convert_pssm_to_csv.py /home/convert_pssm_to_csv.py

COPY run.sh /home/run.sh

CMD /bin/bash /home/run.sh;


