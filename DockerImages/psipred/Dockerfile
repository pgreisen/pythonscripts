FROM ncbi/edirect as edirect
FROM ubuntu:18.04 as blastbuild

# File Author / Maintainer
MAINTAINER Per Greisen <https://github.com/pgreisen>

RUN apt-get -y -m update && apt-get install -y build-essential curl libidn11 libnet-perl perl-doc liblmdb-dev 

RUN apt-get upgrade -y
RUN apt-get autoremove -y
RUN apt-get install tcsh -y
RUN apt-get install -y libgl1-mesa-glx
RUN apt-get install emacs -y
RUN apt-get install unzip zip wget gcc ncbi-blast+ -y


RUN curl -O http://bioinfadmin.cs.ucl.ac.uk/downloads/psipred/psipred.4.02.tar.gz
RUN tar zxf psipred.4.02.tar.gz;
RUN rm psipred.4.02.tar.gz;

RUN cd psipred/src && echo `pwd` && make && \
    make install

RUN wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz && gunzip -v uniref90.fasta.gz
RUN wget http://bioinfadmin.cs.ucl.ac.uk/downloads/pfilt/pfilt1.5.tar.gz && tar zxf pfilt1.5.tar.gz;
RUN cd /pfilt/ && cc -O pfilt.c -lm -o pfilt;
RUN rm pfilt1.5.tar.gz;
RUN /pfilt/pfilt uniref90.fasta > uniref90filt
RUN mkdir uniref90filt_db && mv uniref90filt /uniref90filt_db/uniref90filt_db.fasta
RUN /usr/bin/makeblastdb -dbtype prot -in /uniref90filt_db/uniref90filt_db.fasta -out /uniref90filt_db -parse_seqids -title uniref90filt_db -hash_index
RUN mv uniref90filt_db.* /uniref90filt_db/
RUN rm /uniref90.fasta;