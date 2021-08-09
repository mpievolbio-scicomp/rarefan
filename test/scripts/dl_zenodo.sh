#! /bin/sh

$DATASET=$1

chlororaphis() {
	wget https://zenodo.org/record/5139700/files/chlororaphis.zip?download=1 -O /tmp/zenodo.zip && \
	unzip /tmp/zenodo.zip -d ../data/
}

dokdonia() {
	wget https://zenodo.org/record/5139710/files/dokdonia.zip?download=1 -O /tmp/zenodo.zip && \
	unzip /tmp/zenodo.zip -d ../data/
}

neisseria() {
	wget https://zenodo.org/record/5139705/files/neisseria.zip?download=1 -O /tmp/zenodo.zip && \
	unzip /tmp/zenodo.zip -d ../data/
}
 
all() {
	chlororaphis
	dokdonia
	neisseria
}

"$1"
