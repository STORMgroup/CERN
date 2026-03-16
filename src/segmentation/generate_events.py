#!/usr/bin/env python3
"""
Generate precomputed event values from nanopore signal files using
RawHash or RawHash2 segmentation algorithms.

Reads FAST5, SLOW5/BLOW5, or POD5 files (single files or directories),
processes them through the selected algorithm's event detection pipeline,
and writes a TSV file with read IDs and their event values.

Output format (tab-separated):
    read_id1\tevent1\tevent2\t...\teventn
    read_id2\tevent1\tevent2\t...\teventn

Usage:
    python generate_events.py --method rawhash2 --input /path/to/signals --output events.tsv
    python generate_events.py --method rawhash  --input /path/to/file.fast5 --output events.tsv
"""

import argparse
import sys
from pathlib import Path

import numpy as np

from signal_reader import read_signal_path
from event_detector import (
    detect_events_rawhash,
    detect_events_rawhash2,
    RH1_DEFAULTS,
    RH2_DEFAULTS,
)


def main():
    parser = argparse.ArgumentParser(
        description="Generate precomputed event values from nanopore signal files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input", "-i", required=True,
        help="Path to a signal file (FAST5/POD5/SLOW5/BLOW5) or directory of signal files.",
    )
    parser.add_argument(
        "--output", "-o", required=True,
        help="Output TSV file path.",
    )
    parser.add_argument(
        "--method", "-m", choices=["rawhash", "rawhash2"], default="rawhash2",
        help="Event detection method.",
    )

    # RawHash v1.0 parameters
    rh1 = parser.add_argument_group("RawHash v1.0 parameters")
    rh1.add_argument("--rh1-window1", type=int, default=RH1_DEFAULTS["window_length1"],
                     help="Short window length.")
    rh1.add_argument("--rh1-window2", type=int, default=RH1_DEFAULTS["window_length2"],
                     help="Long window length.")
    rh1.add_argument("--rh1-threshold1", type=float, default=RH1_DEFAULTS["threshold1"],
                     help="Short detector threshold.")
    rh1.add_argument("--rh1-threshold2", type=float, default=RH1_DEFAULTS["threshold2"],
                     help="Long detector threshold.")
    rh1.add_argument("--rh1-peak-height", type=float, default=RH1_DEFAULTS["peak_height"],
                     help="Minimum peak height.")

    # RawHash2 parameters
    rh2 = parser.add_argument_group("RawHash2 parameters")
    rh2.add_argument("--rh2-window1", type=int, default=RH2_DEFAULTS["window_length1"],
                     help="Short window length.")
    rh2.add_argument("--rh2-window2", type=int, default=RH2_DEFAULTS["window_length2"],
                     help="Long window length.")
    rh2.add_argument("--rh2-threshold1", type=float, default=RH2_DEFAULTS["threshold1"],
                     help="Short detector threshold.")
    rh2.add_argument("--rh2-threshold2", type=float, default=RH2_DEFAULTS["threshold2"],
                     help="Long detector threshold.")
    rh2.add_argument("--rh2-peak-height", type=float, default=RH2_DEFAULTS["peak_height"],
                     help="Minimum peak height.")
    rh2.add_argument("--rh2-min-seg-len", type=int, default=RH2_DEFAULTS["min_segment_length"],
                     help="Minimum segment length.")
    rh2.add_argument("--rh2-max-seg-len", type=int, default=RH2_DEFAULTS["max_segment_length"],
                     help="Maximum segment length.")

    args = parser.parse_args()

    # Read signals
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: input path does not exist: {args.input}", file=sys.stderr)
        sys.exit(1)

    print(f"Reading signals from: {args.input}")
    reads = read_signal_path(args.input)
    print(f"  Found {len(reads)} read(s)")

    if not reads:
        print("No reads found. Exiting.", file=sys.stderr)
        sys.exit(1)

    # Process reads
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    n_written = 0
    n_skipped = 0

    with open(output_path, "w") as fout:
        if args.method == "rawhash":
            for read in reads:
                events = detect_events_rawhash(
                    read.signal_pA,
                    window_length1=args.rh1_window1,
                    window_length2=args.rh1_window2,
                    threshold1=args.rh1_threshold1,
                    threshold2=args.rh1_threshold2,
                    peak_height=args.rh1_peak_height,
                )
                if len(events) == 0:
                    n_skipped += 1
                    continue
                event_str = "\t".join(f"{e:.6f}" for e in events)
                fout.write(f"{read.read_id}\t{event_str}\n")
                n_written += 1

        elif args.method == "rawhash2":
            for read in reads:
                events = detect_events_rawhash2(
                    read.signal_pA,
                    window_length1=args.rh2_window1,
                    window_length2=args.rh2_window2,
                    threshold1=args.rh2_threshold1,
                    threshold2=args.rh2_threshold2,
                    peak_height=args.rh2_peak_height,
                    min_segment_length=args.rh2_min_seg_len,
                    max_segment_length=args.rh2_max_seg_len,
                )
                if len(events) == 0:
                    n_skipped += 1
                    continue
                event_str = "\t".join(f"{e:.6f}" for e in events)
                fout.write(f"{read.read_id}\t{event_str}\n")
                n_written += 1

    print(f"Done. Wrote {n_written} read(s) to {args.output} ({n_skipped} skipped)")


if __name__ == "__main__":
    main()
