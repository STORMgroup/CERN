#!/bin/bash
set -e

mkdir -p d1_ecoli_r1041/pod5_files/
cd d1_ecoli_r1041

#Download POD5 from AWS | Unzip; Mv POD5 files into the pod5_files directory. Link: https://github.com/mbhall88/NanoVarBench/blob/main/config/accessions.csv
wget -qO- https://figshare.unimelb.edu.au/ndownloader/files/45408628 | tar xvf -; mv ATCC_25922__202309/*.pod5 pod5_files; rm -rf ATCC_25922__202309
# mv 45408628 ATCC_25922__202309.tar

#Downloading Escherichia coli CFT073, complete genome (https://www.ncbi.nlm.nih.gov/nuccore/AE014075.1/); Unzip; Change name;
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/007/445/GCA_000007445.1_ASM744v1/GCA_000007445.1_ASM744v1_genomic.fna.gz; gunzip GCA_000007445.1_ASM744v1_genomic.fna.gz; mv GCA_000007445.1_ASM744v1_genomic.fna ref.fa

cd ..