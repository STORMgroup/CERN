#!/bin/bash

# This runs the rawhash segmenter on the d9 dataset, and stores the results in "d9_scrappie_events.txt"

OUTDIR="./temp" #Path to the output directory to store all the files to be generated
PREFIX="d1" #A string prefix that you want to attach to the file names to use as an identifier for your current run (e.g., mytestrun)
SIGNALS="../d1_ecoli_r1041/pod5_files" #Path to the directory that contains the fast5 files
REF="../d1_ecoli_r1041/ref.fa" #Path to the reference genome
PORE="../../../extern/uncalled_r1041_model_only_means.txt" #Path to the k-mer model file
PRESETX="sensitive" #Default preset of rawhash for the run (e.g., viral)
THREAD="16" #Number of threads to use
PARAMS="--chunk-size 99999999 -k 9" #(optional -- you can keep it empty) custom parameters to set on top of the default parameters

mkdir temp

/usr/bin/time -vpo "${OUTDIR}/${PREFIX}_rawhash_index_${PRESETX}.time" ../segmentation/rawhash_segment_signal -x ${PRESETX} -t ${THREAD} -p "${PORE}" ${PARAMS} -d "${OUTDIR}/${PREFIX}_rawhash_${PRESETX}.ind" ${REF}
/usr/bin/time -vpo "${OUTDIR}/${PREFIX}_rawhash_map_${PRESETX}.time" ../segmentation/rawhash_segment_signal -x ${PRESETX} -t ${THREAD} ${PARAMS} -o "${OUTDIR}/${PREFIX}_rawhash_${PRESETX}.paf" "${OUTDIR}/${PREFIX}_rawhash_${PRESETX}.ind" ${SIGNALS}

mv events.txt ../d1_ecoli_r1041/d1_scrappieR9_events.txt

sed -i 's/ /\t/g' d1_scrappieR9_events.txt # Replace spaces with tabs, for compatibility with RawHash2

