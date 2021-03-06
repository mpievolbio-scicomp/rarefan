#! /bin/bash

DATASET="$1"

chlororaphis() {
	wget -nv https://zenodo.org/record/5139700/files/chlororaphis.zip?download=1 -O /tmp/zenodo.zip && \
	unzip -q /tmp/zenodo.zip -d ../data/
}

dokdonia() {
	wget -nv https://zenodo.org/record/5139710/files/dokdonia.zip?download=1 -O /tmp/zenodo.zip && \
	unzip -q /tmp/zenodo.zip -d ../data/
}

neisseria() {
	wget -nv https://zenodo.org/record/5139705/files/neisseria.zip?download=1 -O /tmp/zenodo.zip && \
	unzip -q /tmp/zenodo.zip -d ../data/
}
 
all() {
	chlororaphis
	dokdonia
	neisseria
}

"$DATASET"
