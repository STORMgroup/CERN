#!/bin/bash

THREAD=$1

PREFIX="d3_human"
PRESET="fast"
HPC="hpc_off"

bash ../../scripts/run_rawhash2_corrected.sh ${THREAD} ${PREFIX} ${PRESET} ${HPC} scrappieR9

bash ../../scripts/run_rawhash2.sh ${THREAD} ${PREFIX} ${PRESET} ${HPC} scrappieR9

bash ../../scripts/run_rawhash2_corrected.sh ${THREAD} ${PREFIX} ${PRESET} ${HPC} scrappieR10

bash ../../scripts/run_rawhash2.sh ${THREAD} ${PREFIX} ${PRESET} ${HPC} scrappieR10

bash ../../scripts/run_rawhash2_corrected.sh ${THREAD} ${PREFIX} ${PRESET} ${HPC} campolina

bash ../../scripts/run_rawhash2.sh ${THREAD} ${PREFIX} ${PRESET} ${HPC} campolina

