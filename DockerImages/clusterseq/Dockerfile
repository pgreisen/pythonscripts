FROM ubuntu:18.04 

# File Author / Maintainer
MAINTAINER Per Greisen <pgreisen@gmail.com>

RUN apt-get -y -m update && apt-get install -y build-essential curl libidn11 libnet-perl perl-doc liblmdb-dev 

RUN apt-get upgrade -y
RUN apt-get autoremove -y
RUN apt-get install tcsh -y
RUN apt-get install -y libgl1-mesa-glx
RUN apt-get install emacs -y
RUN apt-get install unzip zip wget gcc -y

# Install anaconda
RUN wget https://repo.continuum.io/archive/Anaconda3-2020.02-Linux-x86_64.sh

RUN bash Anaconda3-2020.02-Linux-x86_64.sh -b -p /home/anaconda
RUN rm Anaconda3-2020.02-Linux-x86_64.sh
ENV PATH /home/anaconda/bin:$PATH

COPY requirements.txt /tmp/requirements.txt
RUN pip install --requirement /tmp/requirements.txt

WORKDIR /home/

COPY cluster_sequences_parallel_version.py .

CMD  python cluster_sequences_parallel_version.py -f tot.fasta

