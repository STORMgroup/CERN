#!/bin/bash
set -e

THREAD=$1
PREFIX=$2
PRESET=$3
HPC=$4
SEGMENTER=$5

EXEC="rawhash2"
INDEX="../indexes/${PREFIX}_${HPC}_rawhash2_index_${PRESET}.ind"
SIGNALS="../../data/cern_datasets/CERN_data/${PREFIX}_small/${PREFIX}_small.pod5"

OUTDIR="out_${SEGMENTER}"
EVENTS="../../data/cern_datasets/CERN_data/${PREFIX}_small/${PREFIX}_${SEGMENTER}_events.tsv"
PARAMS="--chunk-size 99999999 --r10 --events-file ${EVENTS}"

if [[ "$HPC" == "hpc_off" ]]; then
    PARAMS="${PARAMS} --sig-diff 0"
fi

mkdir -p "$OUTDIR"

/usr/bin/time -vpo "${OUTDIR}/${PREFIX}_rawhash2_map_${PRESET}.time" \
$EXEC -x ${PRESET} -t ${THREAD} ${PARAMS} \
-o "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.paf" \
"${INDEX}" \
${SIGNALS}