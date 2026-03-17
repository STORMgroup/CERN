#!/bin/bash
set -e

THREAD=32
PORE="../../extern/kmer_models/uncalled_r1041_model_only_means.txt"

OUTDIR="indexes"
mkdir -p "$OUTDIR"

datasets=(d1_ecoli d2_dmelanogaster d3_human)
hpcs=(hpc_on hpc_off)

for d in "${datasets[@]}"; do

  # -------------------------
  # Dataset paths
  # -------------------------
  if [[ "$d" == "d1_ecoli" ]]; then
      PRESET="sensitive"
  elif [[ "$d" == "d2_dmelanogaster" ]]; then
      PRESET="sensitive"
  elif [[ "$d" == "d3_human" ]]; then
      PRESET="fast"
  fi

  REF="../data/cern_datasets/CERN_data/${d}_small/${d}_ref.fa"

  for h in "${hpcs[@]}"; do

    PARAMS="--chunk-size 99999999 --r10"

    if [[ "$h" == "hpc_off" ]]; then
        PARAMS="${PARAMS} --sig-diff -1" # This flag may need to change depending on your version of RawHash2 used
    fi

    PREFIX="${d}_${h}"

    echo "Building RawHash2 index"
    echo "Dataset: $d"
    echo "HPC: $h"

    /usr/bin/time -vpo "${OUTDIR}/${PREFIX}_rawhash2_index_${PRESET}.time" \
    rawhash2 -x ${PRESET} -t ${THREAD} -p "${PORE}" ${PARAMS} \
    -d "${OUTDIR}/${PREFIX}_rawhash2_index_${PRESET}.ind" \
    ${REF}

  done
done