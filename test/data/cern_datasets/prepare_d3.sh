#!/bin/bash
set -e

# This requires Dorado v0.9.0

DATA_DIR="./CERN_data/d3_human_small"
PREFIX="d3_human"

THREAD=$1

# Input / output files
POD5_FILE="${DATA_DIR}/${PREFIX}_small.pod5"
FASTQ_FILE="${DATA_DIR}/${PREFIX}_reads.fastq"
REF_FILE="${DATA_DIR}/${PREFIX}_ref.fa"
PAF_FILE="${DATA_DIR}/${PREFIX}_true_mappings.paf"

# Generate the ground truth mappings:

DORADO_PATH="../../tools/dorado-1.4.0-linux-x64/bin/dorado"

${DORADO_PATH} basecaller sup "$POD5_FILE" --emit-fastq > "$FASTQ_FILE"

minimap2 -x map-ont -t "${THREAD}" -o "$PAF_FILE" "$REF_FILE" "$FASTQ_FILE"

# Generate the segmentations, and correct them:

CERN_PATH="../../../src/run_cern"
SCRAP_PATH="../../../src/segmentation/generate_events.py"
MODEL_DIR="../../../train/models"

# Generate events
python "$SCRAP_PATH" -m rawhash  -i "$POD5_FILE" -o "${DATA_DIR}/${PREFIX}_scrappieR9_events.tsv"

python "$SCRAP_PATH" -m rawhash2 -i "$POD5_FILE" -o "${DATA_DIR}/${PREFIX}_scrappieR10_events.tsv"

# Correct events with CERN
"$CERN_PATH" "$MODEL_DIR/hmm_196_scrapR9.txt"  "${DATA_DIR}/${PREFIX}_scrappieR9_events.tsv"  > "${DATA_DIR}/${PREFIX}_scrappieR9_events_corrected.tsv"

"$CERN_PATH" "$MODEL_DIR/hmm_196_scrapR10.txt" "${DATA_DIR}/${PREFIX}_scrappieR10_events.tsv" > "${DATA_DIR}/${PREFIX}_scrappieR10_events_corrected.tsv"

# Optional: Also create campolina events

# Fill this in:
CAMPOLINA_PATH="../../tools/Campolina"

INFERENCE="${CAMPOLINA_PATH}/inference.py"
SIGNALS="${DATA_DIR}/"
MODEL="${CAMPOLINA_PATH}/R10_model.pth"

python "${INFERENCE}" --pod5_dir "${SIGNALS}" --model_path "${MODEL}" --workers 16 --bs 512

PARQTOE="../../scripts/convert_parquet_to_events.py"
TARGET="${DATA_DIR}/${PREFIX}_campolina_events.tsv"
PARQUET="test_multithread_events.parquet"

python "${PARQTOE}" --parquet "${PARQUET}" --pod5 "${SIGNALS}" --target "${TARGET}"

"$CERN_PATH" "$MODEL_DIR/hmm_196_camp.txt" "$TARGET" > "${DATA_DIR}/${PREFIX}_campolina_events_corrected.tsv"
