FROM conda/miniconda3
MAINTAINER grotec@evolbio.mpg.de

RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install util-linux git make gcc build-essential libgsl-dev gsl-bin andi -y

RUN conda init bash



COPY environment.yml .
RUN conda env create -f environment.yml

SHELL ["conda", "run", "-n", "repinpop", "/bin/bash", "-c"]

