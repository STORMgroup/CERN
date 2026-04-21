#!/bin/bash
set -e

PREFIX=$1
GT=$2
OUTDIR=$3
PAF=$4

# Run pafstats for rawhash2
if ! python ../../scripts/pafstats.py -r "$GT" --annotate "$PAF" \
> "${OUTDIR}/${PREFIX}_rawhash2_ann.paf" \
2> "${OUTDIR}/${PREFIX}_rawhash2.throughput"; then
    echo "Warning: pafstats failed for rawhash2 ($PAF)"
fi