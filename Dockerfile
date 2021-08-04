FROM conda/miniconda3
MAINTAINER grotec@evolbio.mpg.de

# Install system packages.
RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install util-linux git make gcc build-essential libgsl-dev gsl-bin andi -y

# Init conda.
RUN conda init bash

# Copy the environment.yml from the docker build dir into the container's root dir.
COPY environment.yml .

# Create conda env in docker container.
RUN conda env create -f environment.yml
RUN conda activate repinpop

# Build the cmake target
RUN mkdir build
RUN cd build
RUN cmake -DCMAKE_INSTALL_PREFIX ${CONDA_PREFIX} ..
RUN make

# Build java
RUN cd ../REPIN_ecology/REPIN_ecology
RUN gradle build
RUN cd ../../
# SHELL ["conda", "run", "-n", "repinpop", "/bin/bash", "-c"]
