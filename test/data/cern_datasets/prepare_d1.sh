#!/bin/bash
set -e

# This requires Dorado v1.4.0

DATA_DIR="./CERN_data/d1_ecoli_small"
PREFIX="d1_ecoli"

THREAD=$1

# Input / output files
POD5_FILE="${DATA_DIR}/${PREFIX}_small.pod5"
FASTQ_FILE="${DATA_DIR}/${PREFIX}_reads.fastq"
REF_FILE="${DATA_DIR}/${PREFIX}_ref.fa"
PAF_FILE="${DATA_DIR}/${PREFIX}_true_mappings.paf"

# Generate the ground truth mappings:

DORADO_PATH="../../dorado-1.4.0-linux-x64/bin/dorado"

${DORADO_PATH} basecaller sup "$POD5_FILE" --emit-fastq > "$FASTQ_FILE"

minimap2 -x map-ont -t "${THREAD}" -o "$PAF_FILE" "$REF_FILE" "$FASTQ_FILE"

# Generate the segmentations, and correct them:

# Generate events
SEGMENTER="../../../src/segmentation/bin/generate_events"

"$SEGMENTER" -m rawhash  --r10 -i "$POD5_FILE" -o "${DATA_DIR}/${PREFIX}_scrappieR9_events.tsv"

"$SEGMENTER" -m rawhash2 --r10 -i "$POD5_FILE" -o "${DATA_DIR}/${PREFIX}_scrappieR10_events.tsv"


CRANE_PATH="../../../src/crane"
MODEL_DIR="../../../train/models"

# Correct events with CRANE
"$CRANE_PATH" "$MODEL_DIR/hmm_128.tsv" "${DATA_DIR}/${PREFIX}_scrappieR9_events.tsv" 0.3 0.02 -t $THREAD > "${DATA_DIR}/${PREFIX}_scrappieR9_events_corrected.tsv"

"$CRANE_PATH" "$MODEL_DIR/hmm_128.tsv" "${DATA_DIR}/${PREFIX}_scrappieR10_events.tsv" 0.05 0.02 -t $THREAD > "${DATA_DIR}/${PREFIX}_scrappieR10_events_corrected.tsv"

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

"$CRANE_PATH" "$MODEL_DIR/hmm_128.tsv" "$TARGET" 0.035 0.01 -t $THREAD > "${DATA_DIR}/${PREFIX}_campolina_events_corrected.tsv"