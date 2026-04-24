#!/bin/bash
set -e

# This requires Dorado v1.4.0

DATA_DIR="./CERN_data/d1_ecoli_training"
PREFIX="d1_ecoli_training"

THREAD=$1

# Input / output files
FASTQ_FILE="${DATA_DIR}/${PREFIX}_reads.fastq"
REF_FILE="./CERN_data/d1_ecoli_small/d1_ecoli_ref.fa"
PAF_FILE="${DATA_DIR}/${PREFIX}_true_mappings.paf"

# Generate the ground truth mappings:

DORADO_PATH="../../dorado-1.4.0-linux-x64/bin/dorado"

${DORADO_PATH} basecaller sup "$DATA_DIR" --emit-fastq > "$FASTQ_FILE"

minimap2 -x map-ont -t "${THREAD}" -o "$PAF_FILE" "$REF_FILE" "$FASTQ_FILE"

# Generate the segmentations:

# Generate events
SEGMENTER="../../../src/segmentation/bin/generate_events"

"$SEGMENTER" -m rawhash  --r10 -i "$DATA_DIR" -o "${DATA_DIR}/${PREFIX}_scrappieR9_events.tsv"

"$SEGMENTER" -m rawhash2 --r10 -i "$DATA_DIR" -o "${DATA_DIR}/${PREFIX}_scrappieR10_events.tsv"

# Create and correct campolina events

CAMPOLINA_PATH="../../Campolina"

INFERENCE="${CAMPOLINA_PATH}/inference.py"
SIGNALS="${DATA_DIR}/"
MODEL="${CAMPOLINA_PATH}/weights/R10_model.pth"


python "${INFERENCE}" --pod5_dir "${SIGNALS}" --model_path "${MODEL}" --workers 1 --bs 256 --gpu 0

PARQTOE="../../scripts/convert_parquet_to_events.py"
TARGET="${DATA_DIR}/${PREFIX}_campolina_events.tsv"
PARQUET="test_multithread_events.parquet"

python "${PARQTOE}" --parquet "${PARQUET}" --pod5 "${SIGNALS}" --target "${TARGET}"

"$CERN_PATH" "$MODEL_DIR/hmm_196.hmm" "$TARGET" > "${DATA_DIR}/${PREFIX}_campolina_events_corrected.tsv"