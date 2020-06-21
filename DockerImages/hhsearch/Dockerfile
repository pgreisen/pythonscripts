FROM ubuntu:16.04

# Contact me for help!
MAINTAINER p.greisen@enevolv.com

RUN apt-get update && apt-get install -y --no-install-recommends apt-utils
RUN apt-get install -y build-essential;
RUN apt-get install -y emacs wget;
RUN apt-get install -y zlib1g-dev libgl1-mesa-dev;
RUN apt-get install -y awscli;

# Install anaconda
RUN wget https://repo.continuum.io/archive/Anaconda2-2019.10-Linux-x86_64.sh
RUN bash Anaconda2-2019.10-Linux-x86_64.sh -b -p /home/anaconda
RUN rm Anaconda2-2019.10-Linux-x86_64.sh
ENV PATH /home/anaconda/bin:$PATH

COPY requirements.txt /tmp/
RUN pip install --requirement /tmp/requirements.txt

# Copy execution script into container
COPY aws_spinup_and_execute_hhsearch.py /home/aws_spinup_and_execute_hhsearch.py
COPY update_ec2.sh /home/update_ec2.sh

# Files required to run the script
COPY run_hhsearch_aws_nvlv.sh /home/run_hhsearch_aws_nvlv.sh

# insert into configuration file
COPY run.sh /home/run.sh

CMD /bin/bash /home/run.sh;




