#!/bin/bash
set -e

THREAD=$1
MODEL=$2
SEGMENTER=$3

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

PORE="../extern/kmer_models/uncalled_r1041_model_only_means.txt"

mkdir -p "indexes"

REF="../data/cern_datasets/CERN_data/d1_ecoli_small/d1_ecoli_ref.fa"

PARAMS="--chunk-size 99999999 --r10 --sig-diff -1"

if [ -f "indexes/rawhash2_training_index.ind" ]; then
  echo "Index already exists, continuing..."
else
  echo "Creating index..."
  rawhash2 -x sensitive -t ${THREAD} -p "${PORE}" ${PARAMS} \
  -d "./indexes/rawhash2_training_index.ind" \
  ${REF}
fi

# Helper: run the F1 script and extract just the numeric score
get_f1() {
    local pstay=$1
    local pskip=$2
    local winsize=$3
    local result
    result=$(bash correct_and_map_F1.sh "$pstay" "$pskip" "$winsize" "$THREAD" "$MODEL" "$SEGMENTER" 2>/dev/null)
    echo "$result" | grep "F1 Score:" | awk '{print $3}'
}

# Helper: floating point comparison — returns 0 (true) if $1 > $2
fp_gt() {
    awk -v a="$1" -v b="$2" 'BEGIN { exit !(a > b) }'
}

# Helper: add two floats
fp_add() {
    awk -v a="$1" -v b="$2" 'BEGIN { printf "%.4f", a + b }'
}

# Helper: subtract two floats
fp_sub() {
    awk -v a="$1" -v b="$2" 'BEGIN { printf "%.4f", a - b }'
}

echo "========================================"
echo "Phase 1: Searching for best P_STAY and P_SKIP (WINDOWSIZE=20)"
echo "========================================"

WINDOWSIZE=20
PSTAY=0.1
PSKIP=0.1
STEP=0.01

echo "Baseline: PSTAY=$PSTAY, PSKIP=$PSKIP"
BEST_F1=$(get_f1 "$PSTAY" "$PSKIP" "$WINDOWSIZE")
echo "Baseline F1: $BEST_F1"

# Track which direction we just came from to avoid re-testing it
# Format: "pstay_delta,pskip_delta" e.g. "+0.01,0" means we just increased pstay
PREV_PSTAY_DELTA="none"
PREV_PSKIP_DELTA="none"

while true; do
    BEST_PSTAY=$PSTAY
    BEST_PSKIP=$PSKIP
    IMPROVED=false

    # Candidates: (+pstay, 0), (-pstay, 0), (0, +pskip), (0, -pskip)
    # Each entry: "pstay_new pskip_new delta_pstay delta_pskip"
    declare -a CANDIDATES=()

    PSTAY_UP=$(fp_add "$PSTAY" "$STEP")
    PSTAY_DN=$(fp_sub "$PSTAY" "$STEP")
    PSKIP_UP=$(fp_add "$PSKIP" "$STEP")
    PSKIP_DN=$(fp_sub "$PSKIP" "$STEP")

    # Always add all four directions, but we'll skip the one we came from
    # and apply optimization: if we just increased, don't test decreasing (and vice versa)

    # --- PSTAY candidates ---
    # Skip the direction we just came FROM (opposite of last delta)
    ADD_PSTAY_UP=true
    ADD_PSTAY_DN=true

    if [ "$PREV_PSTAY_DELTA" = "up" ]; then
        # We arrived by going up in pstay — came from below, don't go back down
        ADD_PSTAY_DN=false
    elif [ "$PREV_PSTAY_DELTA" = "down" ]; then
        ADD_PSTAY_UP=false
    fi

    # --- PSKIP candidates ---
    ADD_PSKIP_UP=true
    ADD_PSKIP_DN=true

    if [ "$PREV_PSKIP_DELTA" = "up" ]; then
        ADD_PSKIP_DN=false
    elif [ "$PREV_PSKIP_DELTA" = "down" ]; then
        ADD_PSKIP_UP=false
    fi

    # Gather F1s for valid candidates
    BEST_CANDIDATE_F1=$BEST_F1
    BEST_CANDIDATE_PSTAY=$PSTAY
    BEST_CANDIDATE_PSKIP=$PSKIP
    BEST_CANDIDATE_PSTAY_DELTA="none"
    BEST_CANDIDATE_PSKIP_DELTA="none"

    # Test pstay up
    if [ "$ADD_PSTAY_UP" = true ]; then
        echo "  Testing PSTAY=$PSTAY_UP, PSKIP=$PSKIP..."
        F1=$(get_f1 "$PSTAY_UP" "$PSKIP" "$WINDOWSIZE")
        echo "    F1=$F1"
        if fp_gt "$F1" "$BEST_CANDIDATE_F1"; then
            BEST_CANDIDATE_F1=$F1
            BEST_CANDIDATE_PSTAY=$PSTAY_UP
            BEST_CANDIDATE_PSKIP=$PSKIP
            BEST_CANDIDATE_PSTAY_DELTA="up"
            BEST_CANDIDATE_PSKIP_DELTA="none"
        fi
        # Optimization: if increasing pstay helped, skip decreasing pstay
        if fp_gt "$F1" "$BEST_F1"; then
            ADD_PSTAY_DN=false
        fi
    fi

    # Test pstay down
    if [ "$ADD_PSTAY_DN" = true ]; then
        # Guard against going negative
        if awk -v v="$PSTAY_DN" 'BEGIN { exit !(v > 0) }'; then
            echo "  Testing PSTAY=$PSTAY_DN, PSKIP=$PSKIP..."
            F1=$(get_f1 "$PSTAY_DN" "$PSKIP" "$WINDOWSIZE")
            echo "    F1=$F1"
            if fp_gt "$F1" "$BEST_CANDIDATE_F1"; then
                BEST_CANDIDATE_F1=$F1
                BEST_CANDIDATE_PSTAY=$PSTAY_DN
                BEST_CANDIDATE_PSKIP=$PSKIP
                BEST_CANDIDATE_PSTAY_DELTA="down"
                BEST_CANDIDATE_PSKIP_DELTA="none"
            fi
        fi
    fi

    # Test pskip up
    if [ "$ADD_PSKIP_UP" = true ]; then
        echo "  Testing PSTAY=$PSTAY, PSKIP=$PSKIP_UP..."
        F1=$(get_f1 "$PSTAY" "$PSKIP_UP" "$WINDOWSIZE")
        echo "    F1=$F1"
        if fp_gt "$F1" "$BEST_CANDIDATE_F1"; then
            BEST_CANDIDATE_F1=$F1
            BEST_CANDIDATE_PSTAY=$PSTAY
            BEST_CANDIDATE_PSKIP=$PSKIP_UP
            BEST_CANDIDATE_PSTAY_DELTA="none"
            BEST_CANDIDATE_PSKIP_DELTA="up"
        fi
        # Optimization: if increasing pskip helped, skip decreasing pskip
        if fp_gt "$F1" "$BEST_F1"; then
            ADD_PSKIP_DN=false
        fi
    fi

    # Test pskip down
    if [ "$ADD_PSKIP_DN" = true ]; then
        if awk -v v="$PSKIP_DN" 'BEGIN { exit !(v > 0) }'; then
            echo "  Testing PSTAY=$PSTAY, PSKIP=$PSKIP_DN..."
            F1=$(get_f1 "$PSTAY" "$PSKIP_DN" "$WINDOWSIZE")
            echo "    F1=$F1"
            if fp_gt "$F1" "$BEST_CANDIDATE_F1"; then
                BEST_CANDIDATE_F1=$F1
                BEST_CANDIDATE_PSTAY=$PSTAY
                BEST_CANDIDATE_PSKIP=$PSKIP_DN
                BEST_CANDIDATE_PSTAY_DELTA="none"
                BEST_CANDIDATE_PSKIP_DELTA="down"
            fi
        fi
    fi

    # Check if any candidate improved on the current best
    if fp_gt "$BEST_CANDIDATE_F1" "$BEST_F1"; then
        BEST_F1=$BEST_CANDIDATE_F1
        PSTAY=$BEST_CANDIDATE_PSTAY
        PSKIP=$BEST_CANDIDATE_PSKIP
        PREV_PSTAY_DELTA=$BEST_CANDIDATE_PSTAY_DELTA
        PREV_PSKIP_DELTA=$BEST_CANDIDATE_PSKIP_DELTA
        echo "  --> New best: PSTAY=$PSTAY, PSKIP=$PSKIP, F1=$BEST_F1"
    else
        echo "  No improvement found. Converged."
        break
    fi
done

echo "========================================"
echo "Best P_STAY=$PSTAY, P_SKIP=$PSKIP, F1=$BEST_F1 (WINDOWSIZE=$WINDOWSIZE)"
echo "========================================"

echo ""
echo "========================================"
echo "Phase 2: Searching for best WINDOWSIZE (P_STAY=$PSTAY, P_SKIP=$PSKIP)"
echo "========================================"

WIN_STEP=2
PREV_WIN_DELTA="none"

while true; do
    WIN_UP=$((WINDOWSIZE + WIN_STEP))
    WIN_DN=$((WINDOWSIZE - WIN_STEP))

    ADD_WIN_UP=true
    ADD_WIN_DN=true

    if [ "$PREV_WIN_DELTA" = "up" ]; then
        ADD_WIN_DN=false
    elif [ "$PREV_WIN_DELTA" = "down" ]; then
        ADD_WIN_UP=false
    fi

    BEST_WIN=$WINDOWSIZE
    BEST_WIN_F1=$BEST_F1
    BEST_WIN_DELTA="none"

    if [ "$ADD_WIN_UP" = true ]; then
        echo "  Testing WINDOWSIZE=$WIN_UP..."
        F1=$(get_f1 "$PSTAY" "$PSKIP" "$WIN_UP")
        echo "    F1=$F1"
        if fp_gt "$F1" "$BEST_WIN_F1"; then
            BEST_WIN_F1=$F1
            BEST_WIN=$WIN_UP
            BEST_WIN_DELTA="up"
        fi
        # Optimization: if increasing windowsize helped, skip decreasing
        if fp_gt "$F1" "$BEST_F1"; then
            ADD_WIN_DN=false
        fi
    fi

    if [ "$ADD_WIN_DN" = true ] && [ "$WIN_DN" -gt 0 ]; then
        echo "  Testing WINDOWSIZE=$WIN_DN..."
        F1=$(get_f1 "$PSTAY" "$PSKIP" "$WIN_DN")
        echo "    F1=$F1"
        if fp_gt "$F1" "$BEST_WIN_F1"; then
            BEST_WIN_F1=$F1
            BEST_WIN=$WIN_DN
            BEST_WIN_DELTA="down"
        fi
    fi

    if fp_gt "$BEST_WIN_F1" "$BEST_F1"; then
        BEST_F1=$BEST_WIN_F1
        WINDOWSIZE=$BEST_WIN
        PREV_WIN_DELTA=$BEST_WIN_DELTA
        echo "  --> New best: WINDOWSIZE=$WINDOWSIZE, F1=$BEST_F1"
    else
        echo "  No improvement found. Converged."
        break
    fi
done

echo "========================================"
echo "FINAL RESULTS:"
echo "  P_STAY    = $PSTAY"
echo "  P_SKIP    = $PSKIP"
echo "  WINDOWSIZE= $WINDOWSIZE"
echo "  F1 Score  = $BEST_F1"
echo "========================================"