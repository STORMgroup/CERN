#!/bin/bash
set -e

PREFIX="d1_ecoli"
PRESET="sensitive"
HPC="hpc_on"

COMBINED="${PREFIX}_${PRESET}_${HPC}"

OUTDIR="results"
mkdir -p "${OUTDIR}"

TRUE_MAPPINGS="../../data/cern_datasets/CERN_data/${PREFIX}_small/${PREFIX}_true_mappings.paf"

SEGMENTER="scrappieR9"

bash ../../scripts/annotate_paf.sh "${COMBINED}_${SEGMENTER}_corrected" "${TRUE_MAPPINGS}" "${OUTDIR}" "out_${SEGMENTER}_corrected/${PREFIX}_rawhash2_${PRESET}.paf"

bash ../../scripts/annotate_paf.sh "${COMBINED}_${SEGMENTER}" "${TRUE_MAPPINGS}" "${OUTDIR}" "out_${SEGMENTER}/${PREFIX}_rawhash2_${PRESET}.paf"

SEGMENTER="scrappieR10"

bash ../../scripts/annotate_paf.sh "${COMBINED}_${SEGMENTER}_corrected" "${TRUE_MAPPINGS}" "${OUTDIR}" "out_${SEGMENTER}_corrected/${PREFIX}_rawhash2_${PRESET}.paf"

bash ../../scripts/annotate_paf.sh "${COMBINED}_${SEGMENTER}" "${TRUE_MAPPINGS}" "${OUTDIR}" "out_${SEGMENTER}/${PREFIX}_rawhash2_${PRESET}.paf"

SEGMENTER="campolina"

bash ../../scripts/annotate_paf.sh "${COMBINED}_${SEGMENTER}_corrected" "${TRUE_MAPPINGS}" "${OUTDIR}" "out_${SEGMENTER}_corrected/${PREFIX}_rawhash2_${PRESET}.paf"

bash ../../scripts/annotate_paf.sh "${COMBINED}_${SEGMENTER}" "${TRUE_MAPPINGS}" "${OUTDIR}" "out_${SEGMENTER}/${PREFIX}_rawhash2_${PRESET}.paf"