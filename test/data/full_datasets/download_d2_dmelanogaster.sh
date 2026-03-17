#!/bin/bash
set -e

mkdir -p d2_dmelanogaster_r1041/fast5_files/
cd d2_dmelanogaster_r1041

#Download FAST5 from AWS ONT Open Data (D. melanogaster, R10.4.1)
#Source: https://labs.epi2me.io/open-data-dmelanogaster-bkim/
aws s3 cp s3://ont-open-data/contrib/melanogaster_bkim_2023.01/flowcells/D.melanogaster.R1041.400bps/D_melanogaster_1/20221217_1251_MN20261_FAV70669_117da01a/fast5/ ./fast5_files/ --recursive --no-sign-request

#Downloading D. melanogaster Release 6 reference genome (dm6) from UCSC; Unzip; Change name;
wget https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz; gunzip dm6.fa.gz; mv dm6.fa ref.fa

cd ..