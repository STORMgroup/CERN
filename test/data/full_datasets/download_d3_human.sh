#!/bin/bash
set -e

mkdir -p d3_human_hg002_r1041/pod5_files/
cd d3_human_hg002_r1041/pod5_files

#Download POD5 from AWS ONT Open Data (GIAB HG002, R10.4.1)
aws s3 cp s3://ont-open-data/giab_2023.05/flowcells/hg002/20230424_1302_3H_PAO89685_2264ba8c/pod5_pass/ ./ --recursive --no-sign-request

cd ..;

#Downloading CHM13v2 (hs1) Human reference genome;
wget https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz; gunzip hs1.fa.gz; mv hs1.fa ref.fa

cd ..