#!/bin/bash
set -e

THREAD=$1
MODEL=$2
SEGMENTER=$3

# Derive a unique tag for this run so parallel runs don't collide
MODEL_TAG=$(basename "$MODEL" | sed 's/\.[^.]*$//')
RUN_ID="${MODEL_TAG}_${SEGMENTER}"

echo "========================================"
echo "Parameter search starting"
echo "  Model    : $MODEL  ($MODEL_TAG)"
echo "  Segmenter: $SEGMENTER"
echo "========================================"
echo ""

cd ../../test/data/cern_datasets

if [ -d "CERN_data" ]; then
    echo "CERN_data directory exists, continuing..."
else
    echo "CERN_data not found, running download script..."
    bash download_cern_data.sh
fi

if [ -f "CERN_data/d1_ecoli_training/d1_ecoli_training_true_mappings.paf" ]; then
  echo "Ground truth exists, continuing..."
else
  echo "Please run prepare_training.sh with GPU access"
  exit 1
fi

cd ../../../train/tuning_parameters

DATA_DIR="../../test/data/cern_datasets/CERN_data/d1_ecoli_training/"

PORE="../../extern/kmer_models/uncalled_r1041_model_only_means.txt"

mkdir -p "indexes"

REF="../../test/data/cern_datasets/CERN_data/d1_ecoli_small/d1_ecoli_ref.fa"

PARAMS="--chunk-size 99999999 --r10 --sig-diff -1"

if [ -f "indexes/rawhash2_training_index.ind" ]; then
  echo "Index already exists, continuing..."
else
  echo "Creating index..."
  rawhash2 -x sensitive -t ${THREAD} -p "${PORE}" ${PARAMS} \
  -d "./indexes/rawhash2_training_index.ind" \
  ${REF}
fi

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

get_f1() {
    local pstay=$1 pskip=$2 winsize=$3
    local result
    result=$(bash correct_and_map_F1.sh "$pstay" "$pskip" "$winsize" "$THREAD" "$MODEL" "$SEGMENTER" "$RUN_ID" 2>/dev/null)
    echo "$result" | grep "F1 Score:" | awk '{print $3}'
}

fp_gt() { awk -v a="$1" -v b="$2" 'BEGIN { exit !(a > b) }'; }
fp_add() { awk -v a="$1" -v b="$2" 'BEGIN { printf "%.6f", a + b }'; }
fp_sub() { awk -v a="$1" -v b="$2" 'BEGIN { printf "%.6f", a - b }'; }
fp_pos() { awk -v a="$1" 'BEGIN { exit !(a >= 0) }'; }  # returns true if >= 0

# ---------------------------------------------------------------------------
# Hill-climb for pstay/pskip at a fixed windowsize and starting point,
# with a given step size. Updates globals PSTAY, PSKIP, BEST_F1.
# ---------------------------------------------------------------------------
hillclimb_pstay_pskip() {
    local step=$1

    local PREV_PSTAY_DELTA="none"
    local PREV_PSKIP_DELTA="none"

    echo "  [Hill-climb] Starting from PSTAY=$PSTAY PSKIP=$PSKIP F1=$BEST_F1 step=$step"

    while true; do
        local PSTAY_UP PSTAY_DN PSKIP_UP PSKIP_DN
        PSTAY_UP=$(fp_add "$PSTAY" "$step")
        PSTAY_DN=$(fp_sub "$PSTAY" "$step")
        PSKIP_UP=$(fp_add "$PSKIP" "$step")
        PSKIP_DN=$(fp_sub "$PSKIP" "$step")

        local ADD_PSTAY_UP=true ADD_PSTAY_DN=true
        local ADD_PSKIP_UP=true ADD_PSKIP_DN=true

        [ "$PREV_PSTAY_DELTA" = "up"   ] && ADD_PSTAY_DN=false
        [ "$PREV_PSTAY_DELTA" = "down" ] && ADD_PSTAY_UP=false
        [ "$PREV_PSKIP_DELTA" = "up"   ] && ADD_PSKIP_DN=false
        [ "$PREV_PSKIP_DELTA" = "down" ] && ADD_PSKIP_UP=false

        local BEST_CANDIDATE_F1=$BEST_F1
        local BEST_CANDIDATE_PSTAY=$PSTAY
        local BEST_CANDIDATE_PSKIP=$PSKIP
        local BEST_CANDIDATE_PSTAY_DELTA="none"
        local BEST_CANDIDATE_PSKIP_DELTA="none"

        # Test pstay up
        if [ "$ADD_PSTAY_UP" = true ]; then
            echo "    Testing PSTAY=$PSTAY_UP PSKIP=$PSKIP..."
            local F1; F1=$(get_f1 "$PSTAY_UP" "$PSKIP" "$WINDOWSIZE")
            echo "      F1=$F1"
            if fp_gt "$F1" "$BEST_CANDIDATE_F1"; then
                BEST_CANDIDATE_F1=$F1
                BEST_CANDIDATE_PSTAY=$PSTAY_UP
                BEST_CANDIDATE_PSKIP=$PSKIP
                BEST_CANDIDATE_PSTAY_DELTA="up"
                BEST_CANDIDATE_PSKIP_DELTA="none"
            fi
            fp_gt "$F1" "$BEST_F1" && ADD_PSTAY_DN=false
        fi

        # Test pstay down
        if [ "$ADD_PSTAY_DN" = true ] && fp_pos "$PSTAY_DN"; then
            echo "    Testing PSTAY=$PSTAY_DN PSKIP=$PSKIP..."
            local F1; F1=$(get_f1 "$PSTAY_DN" "$PSKIP" "$WINDOWSIZE")
            echo "      F1=$F1"
            if fp_gt "$F1" "$BEST_CANDIDATE_F1"; then
                BEST_CANDIDATE_F1=$F1
                BEST_CANDIDATE_PSTAY=$PSTAY_DN
                BEST_CANDIDATE_PSKIP=$PSKIP
                BEST_CANDIDATE_PSTAY_DELTA="down"
                BEST_CANDIDATE_PSKIP_DELTA="none"
            fi
        fi

        # Test pskip up
        if [ "$ADD_PSKIP_UP" = true ]; then
            echo "    Testing PSTAY=$PSTAY PSKIP=$PSKIP_UP..."
            local F1; F1=$(get_f1 "$PSTAY" "$PSKIP_UP" "$WINDOWSIZE")
            echo "      F1=$F1"
            if fp_gt "$F1" "$BEST_CANDIDATE_F1"; then
                BEST_CANDIDATE_F1=$F1
                BEST_CANDIDATE_PSTAY=$PSTAY
                BEST_CANDIDATE_PSKIP=$PSKIP_UP
                BEST_CANDIDATE_PSTAY_DELTA="none"
                BEST_CANDIDATE_PSKIP_DELTA="up"
            fi
            fp_gt "$F1" "$BEST_F1" && ADD_PSKIP_DN=false
        fi

        # Test pskip down
        if [ "$ADD_PSKIP_DN" = true ] && fp_pos "$PSKIP_DN"; then
            echo "    Testing PSTAY=$PSTAY PSKIP=$PSKIP_DN..."
            local F1; F1=$(get_f1 "$PSTAY" "$PSKIP_DN" "$WINDOWSIZE")
            echo "      F1=$F1"
            if fp_gt "$F1" "$BEST_CANDIDATE_F1"; then
                BEST_CANDIDATE_F1=$F1
                BEST_CANDIDATE_PSTAY=$PSTAY
                BEST_CANDIDATE_PSKIP=$PSKIP_DN
                BEST_CANDIDATE_PSTAY_DELTA="none"
                BEST_CANDIDATE_PSKIP_DELTA="down"
            fi
        fi

        if fp_gt "$BEST_CANDIDATE_F1" "$BEST_F1"; then
            BEST_F1=$BEST_CANDIDATE_F1
            PSTAY=$BEST_CANDIDATE_PSTAY
            PSKIP=$BEST_CANDIDATE_PSKIP
            PREV_PSTAY_DELTA=$BEST_CANDIDATE_PSTAY_DELTA
            PREV_PSKIP_DELTA=$BEST_CANDIDATE_PSKIP_DELTA
            echo "    --> New best: PSTAY=$PSTAY PSKIP=$PSKIP F1=$BEST_F1"
        else
            echo "    Converged at step=$step"
            break
        fi
    done
}

# ---------------------------------------------------------------------------
# Phase 1: Grid search over {0, 0.1, 0.2, 0.3} x {0, 0.1, 0.2, 0.3}
# ---------------------------------------------------------------------------
echo "========================================"
echo "Phase 1: Grid search (WINDOWSIZE=20)"
echo "========================================"

WINDOWSIZE=20
GRID_VALS="0.0 0.1 0.2 0.3"

BEST_F1="-1"
PSTAY="0.1"
PSKIP="0.1"

for gs in $GRID_VALS; do
    for gk in $GRID_VALS; do
        echo "  Grid point PSTAY=$gs PSKIP=$gk..."
        F1=$(get_f1 "$gs" "$gk" "$WINDOWSIZE")
        echo "    F1=$F1"
        if fp_gt "$F1" "$BEST_F1"; then
            BEST_F1=$F1
            PSTAY=$gs
            PSKIP=$gk
            echo "    --> New grid best!"
        fi
    done
done

echo ""
echo "Grid search best: PSTAY=$PSTAY PSKIP=$PSKIP F1=$BEST_F1"

# ---------------------------------------------------------------------------
# Phase 2: Multi-scale hill-climb for pstay/pskip
# ---------------------------------------------------------------------------
echo ""
echo "========================================"
echo "Phase 2: Multi-scale hill-climb for P_STAY and P_SKIP (WINDOWSIZE=20)"
echo "========================================"

for STEP in 0.02 0.01 0.005; do
    hillclimb_pstay_pskip "$STEP"
done

echo ""
echo "Best P_STAY=$PSTAY P_SKIP=$PSKIP F1=$BEST_F1 (WINDOWSIZE=$WINDOWSIZE)"

# ---------------------------------------------------------------------------
# Phase 3: Hill-climb for WINDOWSIZE
# ---------------------------------------------------------------------------
echo ""
echo "========================================"
echo "Phase 3: Hill-climb for WINDOWSIZE (P_STAY=$PSTAY P_SKIP=$PSKIP)"
echo "========================================"

WIN_STEP=2
PREV_WIN_DELTA="none"

while true; do
    WIN_UP=$((WINDOWSIZE + WIN_STEP))
    WIN_DN=$((WINDOWSIZE - WIN_STEP))

    ADD_WIN_UP=true
    ADD_WIN_DN=true
    [ "$PREV_WIN_DELTA" = "up"   ] && ADD_WIN_DN=false
    [ "$PREV_WIN_DELTA" = "down" ] && ADD_WIN_UP=false

    BEST_WIN=$WINDOWSIZE
    BEST_WIN_F1=$BEST_F1
    BEST_WIN_DELTA="none"

    if [ "$ADD_WIN_UP" = true ]; then
        echo "  Testing WINDOWSIZE=$WIN_UP..."
        F1=$(get_f1 "$PSTAY" "$PSKIP" "$WIN_UP")
        echo "    F1=$F1"
        if fp_gt "$F1" "$BEST_WIN_F1"; then
            BEST_WIN_F1=$F1; BEST_WIN=$WIN_UP; BEST_WIN_DELTA="up"
        fi
        fp_gt "$F1" "$BEST_F1" && ADD_WIN_DN=false
    fi

    if [ "$ADD_WIN_DN" = true ] && [ "$WIN_DN" -gt 0 ]; then
        echo "  Testing WINDOWSIZE=$WIN_DN..."
        F1=$(get_f1 "$PSTAY" "$PSKIP" "$WIN_DN")
        echo "    F1=$F1"
        if fp_gt "$F1" "$BEST_WIN_F1"; then
            BEST_WIN_F1=$F1; BEST_WIN=$WIN_DN; BEST_WIN_DELTA="down"
        fi
    fi

    if fp_gt "$BEST_WIN_F1" "$BEST_F1"; then
        BEST_F1=$BEST_WIN_F1
        WINDOWSIZE=$BEST_WIN
        PREV_WIN_DELTA=$BEST_WIN_DELTA
        echo "  --> New best: WINDOWSIZE=$WINDOWSIZE F1=$BEST_F1"
    else
        echo "  Converged."
        break
    fi
done

# ---------------------------------------------------------------------------
# Final results
# ---------------------------------------------------------------------------
echo ""
echo "========================================"
echo "FINAL RESULTS:"
echo "  Model     = $MODEL"
echo "  Segmenter = $SEGMENTER"
echo "  P_STAY    = $PSTAY"
echo "  P_SKIP    = $PSKIP"
echo "  WINDOWSIZE= $WINDOWSIZE"
echo "  F1 Score  = $BEST_F1"
echo "========================================"