#!/usr/bin/env python3
"""
pafstats.py — Evaluation script for RawHash2 read mapping accuracy.

Compares RawHash2 PAF output against a minimap2 ground-truth PAF.
Uses overlap-based TP/FP/FN/TN classification (same logic as UNCALLED pafstats).

Usage (compatible with 'uncalled pafstats' interface):
    python3 pafstats.py -r <ground_truth.paf> --annotate <rawhash2.paf> \
        > annotated.paf 2> throughput.txt

Output:
    stdout: Annotated PAF (each line gets rf:Z:tp/fp/fn/tn/na appended)
    stderr: Summary stats, confusion matrix, throughput (BP/sec)

The annotated PAF is compatible with analyze_paf.py (expects 20 cols for mapped,
15 cols for unmapped reads).

Overlap-based classification logic from UNCALLED:
    https://github.com/skovaka/UNCALLED/blob/master/uncalled/pafstats.py
"""

import sys
import argparse
import numpy as np


class PafEntry:
    """Represents a single PAF alignment record."""

    def __init__(self, line):
        tabs = line.rstrip().split('\t')
        self.qr_name = tabs[0]
        self.qr_len = int(tabs[1])
        self.is_mapped = tabs[4] != '*'

        if self.is_mapped:
            self.qr_st = int(tabs[2])
            self.qr_en = int(tabs[3])
            self.is_fwd = tabs[4] == '+'
            self.rf_name = tabs[5]
            self.rf_len = int(tabs[6])
            self.rf_st = int(tabs[7])
            self.rf_en = int(tabs[8])
        else:
            self.qr_st = 1
            self.qr_en = self.qr_len
            self.is_fwd = None
            self.rf_name = None
            self.rf_len = None
            self.rf_st = None
            self.rf_en = None

        self.raw_line = line.rstrip()

    def ext_ref(self, ext=1.0):
        """Extend reference coordinates by a factor based on unmapped query ends."""
        st_shift = int(self.qr_st * ext)
        en_shift = int((self.qr_len - self.qr_en) * ext)

        if self.is_fwd:
            return (max(1, self.rf_st - st_shift),
                    min(self.rf_len, self.rf_en + en_shift))
        else:
            return (max(1, self.rf_st - en_shift),
                    min(self.rf_len, self.rf_en + st_shift))

    def overlaps(self, other, ext=0.0):
        """Check if two mapped entries overlap on the reference (with extension)."""
        if not (self.is_mapped and other.is_mapped):
            return False
        if self.rf_name != other.rf_name:
            return False
        st1, en1 = self.ext_ref(ext)
        st2, en2 = other.ext_ref(ext)
        return max(st1, st2) <= min(en1, en2)


def parse_paf(filepath):
    """Parse a PAF file, yielding PafEntry objects. Skips comment/log lines."""
    with open(filepath) as f:
        for line in f:
            if line.startswith('#') or line.startswith('['):
                continue
            tabs = line.rstrip().split('\t')
            if len(tabs) < 12:
                # Skip malformed/truncated lines
                continue
            yield PafEntry(line)


def classify_reads(query_entries, ref_entries, ext=1.5):
    """
    Classify query reads against reference using overlap-based comparison.

    For each query read:
      - If mapped in query:
        - If reference has no mapping for this read → 'na' (false positive, unmapped in ref)
        - If overlaps reference mapping (with ext) → 'tp'
        - Otherwise → 'fp'
      - If unmapped in query:
        - If reference has no mapping → 'tn'
        - If reference has a mapping → 'fn'

    Returns list of (PafEntry, label) tuples.
    """
    # Build reference lookup: read_name → list of PafEntry
    ref_locs = {}
    for r in ref_entries:
        ref_locs.setdefault(r.qr_name, []).append(r)

    results = []

    for q in query_entries:
        refs = ref_locs.get(q.qr_name, None)

        if q.is_mapped:
            if refs is None or not refs[0].is_mapped:
                # Mapped in query but not in reference
                results.append((q, 'na'))
                continue

            match = False
            for r in refs:
                if q.overlaps(r, ext):
                    match = True
                    break

            if match:
                results.append((q, 'tp'))
            else:
                results.append((q, 'fp'))
        else:
            if refs is None or not refs[0].is_mapped:
                results.append((q, 'tn'))
            else:
                results.append((q, 'fn'))

    return results


def main():
    parser = argparse.ArgumentParser(
        description='Compare RawHash2 PAF against minimap2 ground truth.')
    parser.add_argument('infile', type=str,
                        help='RawHash2 PAF file')
    parser.add_argument('-r', '--ref-paf', required=True, type=str,
                        help='Reference (minimap2) PAF file')
    parser.add_argument('-a', '--annotate', action='store_true',
                        help='Output annotated PAF to stdout (with rf:Z: tag)')
    args = parser.parse_args()

    statsout = sys.stderr if args.annotate else sys.stdout

    # Parse both PAF files
    query_entries = list(parse_paf(args.infile))
    ref_entries = list(parse_paf(args.ref_paf))

    num_mapped = sum(1 for p in query_entries if p.is_mapped)
    statsout.write("Summary: %d reads, %d mapped (%.2f%%)\n\n" %
                   (len(query_entries), num_mapped,
                    100 * num_mapped / len(query_entries) if query_entries else 0))

    # Classify reads
    results = classify_reads(query_entries, ref_entries)

    # Count classifications
    counts = {'tp': 0, 'tn': 0, 'fp': 0, 'fn': 0, 'na': 0}
    for _, label in results:
        counts[label] += 1

    n = len(results)
    ntp, ntn, nfp, nfn, nna = counts['tp'], counts['tn'], counts['fp'], counts['fn'], counts['na']

    statsout.write("Comparing to reference PAF\n")
    statsout.write("     P     N\n")
    statsout.write("T %6.2f %5.2f\n" % (100 * ntp / n, 100 * ntn / n))
    statsout.write("F %6.2f %5.2f\n" % (100 * (nfp) / n, 100 * nfn / n))
    statsout.write("NA: %.2f\n\n" % (100 * nna / n))

    # Output annotated PAF
    if args.annotate:
        for entry, label in results:
            sys.stdout.write("%s\trf:Z:%s\n" % (entry.raw_line, label))

    # Throughput stats (from mt:f: tag)
    map_ms = []
    map_bp = []
    for entry, _ in results:
        if entry.is_mapped:
            # Extract mt:f: tag from raw line
            for field in entry.raw_line.split('\t')[12:]:
                if field.startswith('mt:f:'):
                    mt = float(field.split(':')[2])
                    map_ms.append(mt)
                    map_bp.append(entry.qr_en)
                    break

    if map_ms:
        map_ms = np.array(map_ms)
        map_bp = np.array(map_bp)
        map_bpps = 1000 * map_bp / map_ms

        statsout.write("Speed            Mean    Median\n")
        statsout.write("BP per sec: %9.2f %9.2f\n" % (np.mean(map_bpps), np.median(map_bpps)))
        statsout.write("BP mapped:  %9.2f %9.2f\n" % (np.mean(map_bp), np.median(map_bp)))
        statsout.write("MS to map:  %9.2f %9.2f\n" % (np.mean(map_ms), np.median(map_ms)))


if __name__ == '__main__':
    main()