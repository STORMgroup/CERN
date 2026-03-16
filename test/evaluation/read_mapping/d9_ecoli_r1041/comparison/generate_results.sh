#!/bin/bash

uncalled pafstats -r ../true_mappings.paf --annotate ../rawhash2/d1_ecoli_r1041_rawhash2_sensitive.paf > d1_ecoli_r1041_rawhash2_sensitive_ann.paf 2> d1_ecoli_r1041_rawhash2_sensitive.throughput
uncalled pafstats -r ../true_mappings.paf --annotate ../rawhash2corrected/d1_ecoli_r1041_rawhash2_sensitive.paf > d1_ecoli_r1041_rawhash2_corrected_sensitive_ann.paf 2> d1_ecoli_r1041_rawhash2_corrected_sensitive.throughput

python ../../../../scripts/analyze_paf.py d1_ecoli_r1041_rawhash2_sensitive_ann.paf > d1_ecoli_r1041_rawhash2_results.txt
python ../../../../scripts/analyze_paf.py d1_ecoli_r1041_rawhash2_sensitive_ann.paf > d1_ecoli_r1041_rawhash2_corrected_results.txt

echo "Default rawhash2:"
cat d1_ecoli_r1041_rawhash2_results.txt

echo "rawhash2 with corrected events:"
cat d1_ecoli_r1041_rawhash2_corrected_results.txt
