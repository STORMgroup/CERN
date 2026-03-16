#!/bin/bash

THREAD=$1

#d9_ecoli_r1041
OUTDIR="./rawhash2/"
FAST5="../../../data/d1_ecoli_r1041/pod5_files/"
REF="../../../data/d1_ecoli_r1041/ref.fa"
PORE="../../../../extern/kmer_models/uncalled_r1041_model_only_means.txt"
PRESET="sensitive"
mkdir -p ${OUTDIR}
PARAMS="--r10"

#The following is the run using default parameters:
PREFIX="d1_ecoli_r1041"
bash ../../../scripts/run_rawhash2.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${PORE} ${PRESET} ${THREAD} "${PARAMS}" > "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.out" 2> "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.err"
