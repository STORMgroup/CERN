#!/bin/bash
set -e

OUTFILE="summary.txt"
> "$OUTFILE"   # clear previous output

find . -type f -name "*.throughput" | sort | while read -r file; do
    echo "File: $file" >> "$OUTFILE"

    if grep -q "Traceback" "$file"; then
        echo "out of memory" >> "$OUTFILE"
        echo "" >> "$OUTFILE"
        continue
    fi

    awk '
    /Summary:/ {
        reads=$2
        mapped=$4
        gsub(",", "", reads)
        gsub(",", "", mapped)
        summary=$0
    }

    /^T/ {
        TP=$2
        TN=$3
    }

    /^F/ {
        FP=$2
        FN=$3
    }

    /^NA:/ {
        NA=$2
    }

    /BP per sec:/ {
        bps_mean=$4
        bps_median=$5
    }

    /BP mapped:/ {
        bpm_mean=$3
        bpm_median=$4
    }

    /MS to map:/ {
        ms_mean=$4
        ms_median=$5
    }

    END {
        precision = (TP + FP > 0) ? TP / (TP + FP) : 0
        recall = (TP + FN > 0) ? TP / (TP + FN) : 0
        f1 = (precision + recall > 0) ? 2 * precision * recall / (precision + recall) : 0
        accuracy = (TP + TN) / (TP + TN + FP + FN)

        printf "%s\n", summary
        printf "TP: %.2f  FP: %.2f  FN: %.2f  TN: %.2f\n", TP, FP, FN, TN
        printf "Precision: %.4f\n", precision
        printf "Recall: %.4f\n", recall
        printf "F1 score: %.4f\n", f1
        printf "Accuracy: %.4f\n", accuracy
        printf "BP per sec (mean/median): %s / %s\n", bps_mean, bps_median
        printf "BP mapped (mean/median): %s / %s\n", bpm_mean, bpm_median
        printf "MS to map (mean/median): %s / %s\n", ms_mean, ms_median
        print ""
    }
    ' "$file" >> "$OUTFILE"

done

echo "Results written to $OUTFILE"
