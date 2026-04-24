#!/bin/bash
set -e

cd ../../test/data/cern_datasets

if [ -d "CERN_data" ]; then
    echo "CERN_data directory exists, continuing..."
else
    echo "CERN_data not found, running download script..."
    bash download_cern_data.sh
fi

bash prepare_training.sh

cd ../../../train/tuning_params

DATA_DIR="../../test/data/cern_datasets/CERN_data/d1_ecoli_training/"





THREAD=32
PORE="../../extern/kmer_models/uncalled_r1041_model_only_means.txt"

OUTDIR="indexes"
mkdir -p "$OUTDIR"


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