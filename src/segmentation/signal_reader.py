"""
Signal reader for nanopore signal files (FAST5, SLOW5/BLOW5, POD5).

Reads raw signals, calibrates them to picoamperes (pA), and applies
the same 30-200 pA hardware filter used by RawHash/RawHash2.
"""

import os
from pathlib import Path
from dataclasses import dataclass

import numpy as np


@dataclass
class Read:
    """A single nanopore read with calibrated signal."""
    read_id: str
    signal_pA: np.ndarray  # float32, pA values after 30-200 filter


def _filter_pA(signal: np.ndarray) -> np.ndarray:
    """Apply the 30-200 pA hardware filter (matches RawHash/RawHash2 rsig.c)."""
    mask = (signal > 30.0) & (signal < 200.0)
    return signal[mask].astype(np.float32)


# ── FAST5 reader ─────────────────────────────────────────────────────────────

def _read_fast5(filepath: str) -> list[Read]:
    """Read all reads from a FAST5 file."""
    import h5py

    reads = []
    with h5py.File(filepath, "r") as f:
        # Handle both single-read and multi-read FAST5
        if "Raw" in f:
            # Single-read FAST5
            read_groups = [f]
        else:
            # Multi-read FAST5: each read is under /read_XXXX
            read_groups = [f[k] for k in f.keys() if k.startswith("read_")]

        for rg in read_groups:
            # Find the read group
            if "Raw" in rg:
                raw_group = rg["Raw"]
            else:
                continue

            # Get read_id
            if "Reads" in raw_group:
                # Multi-read format: Raw/Reads/Read_XXXX
                for read_key in raw_group["Reads"]:
                    read_data = raw_group["Reads"][read_key]
                    read_id = read_data.attrs.get("read_id", read_key)
                    if isinstance(read_id, bytes):
                        read_id = read_id.decode()
                    raw_signal = read_data["Signal"][:]

                    # Get channel info for calibration
                    channel = rg["UniqueGlobalKey/channel_id"]
                    digitisation = float(channel.attrs["digitisation"])
                    offset = float(channel.attrs["offset"])
                    ran = float(channel.attrs["range"])
                    scale = ran / digitisation

                    pA = (raw_signal.astype(np.float32) + offset) * scale
                    reads.append(Read(read_id=read_id, signal_pA=_filter_pA(pA)))
            else:
                # Single-read: Raw/Signal directly, or read attrs
                read_id = raw_group.attrs.get("read_id", "unknown")
                if isinstance(read_id, bytes):
                    read_id = read_id.decode()

                if "Signal" in raw_group:
                    raw_signal = raw_group["Signal"][:]
                else:
                    continue

                channel = rg["UniqueGlobalKey/channel_id"]
                digitisation = float(channel.attrs["digitisation"])
                offset = float(channel.attrs["offset"])
                ran = float(channel.attrs["range"])
                scale = ran / digitisation

                pA = (raw_signal.astype(np.float32) + offset) * scale
                reads.append(Read(read_id=read_id, signal_pA=_filter_pA(pA)))

    return reads


# ── POD5 reader ──────────────────────────────────────────────────────────────

def _read_pod5(filepath: str) -> list[Read]:
    """Read all reads from a POD5 file."""
    import pod5

    reads = []
    with pod5.Reader(filepath) as reader:
        for read in reader.reads():
            read_id = str(read.read_id)
            raw_signal = read.signal
            cal = read.calibration
            pA = (raw_signal.astype(np.float32) + cal.offset) * cal.scale
            reads.append(Read(read_id=read_id, signal_pA=_filter_pA(pA)))

    return reads


# ── SLOW5/BLOW5 reader ───────────────────────────────────────────────────────

def _read_slow5(filepath: str) -> list[Read]:
    """Read all reads from a SLOW5 or BLOW5 file."""
    import pyslow5

    reads = []
    s5 = pyslow5.Open(filepath, "r")
    read_list = s5.seq_reads()
    for rec in read_list:
        read_id = rec["read_id"]
        raw_signal = rec["signal"]
        digitisation = float(rec["digitisation"])
        offset = float(rec["offset"])
        ran = float(rec["range"])
        scale = ran / digitisation
        pA = (raw_signal.astype(np.float32) + offset) * scale
        reads.append(Read(read_id=read_id, signal_pA=_filter_pA(pA)))

    return reads


# ── Unified interface ────────────────────────────────────────────────────────

_READERS = {
    ".fast5": _read_fast5,
    ".pod5": _read_pod5,
    ".slow5": _read_slow5,
    ".blow5": _read_slow5,
}


def read_signal_file(filepath: str) -> list[Read]:
    """Read all reads from a signal file (FAST5, POD5, SLOW5, BLOW5)."""
    ext = Path(filepath).suffix.lower()
    reader = _READERS.get(ext)
    if reader is None:
        raise ValueError(f"Unsupported file format: {ext} (supported: {list(_READERS.keys())})")
    return reader(filepath)


def read_signal_path(path: str) -> list[Read]:
    """Read signals from a single file or all signal files in a directory."""
    p = Path(path)
    if p.is_file():
        return read_signal_file(str(p))

    if p.is_dir():
        reads = []
        for f in sorted(p.iterdir()):
            if f.suffix.lower() in _READERS:
                try:
                    reads.extend(read_signal_file(str(f)))
                except Exception as e:
                    print(f"Warning: failed to read {f}: {e}")
        return reads

    raise FileNotFoundError(f"Path does not exist: {path}")
