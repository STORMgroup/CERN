#!/bin/bash

THREAD=$1

bash ../../../scripts/run_minimap2.sh . ../../../data/d1_ecoli_r1041/reads.fastq ../../../data/d1_ecoli_r1041/ref.fa ${THREAD}
