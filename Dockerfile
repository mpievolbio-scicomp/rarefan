FROM conda/miniconda3
MAINTAINER computing@evolbio.mpg.de

# Install system packages.
RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install linux-libc-dev
RUN apt-get install util-linux git make gcc build-essential libgsl-dev gsl-bin andi wget zip unzip -y

# Init conda.
RUN conda init bash

# Copy the environment.yml from the docker build dir into the container's root dir.
COPY environment.yml .
COPY CI .

# Create conda env in docker container.
RUN conda env create -f environment.yml
