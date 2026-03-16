"""
FFI wrapper for RawHash and RawHash2 event detection shared libraries.

Uses ctypes to call the C implementations, preserving the exact algorithms
from the original RawHash (v1.0) and RawHash2 codebases.
"""

import ctypes
from pathlib import Path

import numpy as np

_LIB_DIR = Path(__file__).parent / "lib"


def _load_lib(name: str) -> ctypes.CDLL:
    """Load a shared library from the lib/ directory."""
    so_path = _LIB_DIR / name
    if not so_path.exists():
        raise FileNotFoundError(
            f"Shared library not found: {so_path}\n"
            f"Run 'make' in {_LIB_DIR.parent} to compile."
        )
    return ctypes.CDLL(str(so_path))


# ── Lazy-loaded library handles ──────────────────────────────────────────────

_RH1_LIB = None
_RH2_LIB = None


def _get_rh1_lib() -> ctypes.CDLL:
    global _RH1_LIB
    if _RH1_LIB is not None:
        return _RH1_LIB
    _RH1_LIB = _load_lib("rawhash_events.so")

    _RH1_LIB.rh1_detect_events.restype = ctypes.POINTER(ctypes.c_float)
    _RH1_LIB.rh1_detect_events.argtypes = [
        ctypes.c_uint32,  # s_len
        ctypes.POINTER(ctypes.c_float),  # sig
        ctypes.POINTER(ctypes.c_uint32),  # n_events (out)
        ctypes.c_uint32,  # window_length1
        ctypes.c_uint32,  # window_length2
        ctypes.c_float,  # threshold1
        ctypes.c_float,  # threshold2
        ctypes.c_float,  # peak_height
    ]
    _RH1_LIB.rh1_free.argtypes = [ctypes.c_void_p]
    _RH1_LIB.rh1_free.restype = None

    return _RH1_LIB


def _get_rh2_lib() -> ctypes.CDLL:
    global _RH2_LIB
    if _RH2_LIB is not None:
        return _RH2_LIB
    _RH2_LIB = _load_lib("rawhash2_events.so")

    _RH2_LIB.rh2_detect_events.restype = ctypes.POINTER(ctypes.c_float)
    _RH2_LIB.rh2_detect_events.argtypes = [
        ctypes.c_uint32,  # s_len
        ctypes.POINTER(ctypes.c_float),  # sig
        ctypes.POINTER(ctypes.c_uint32),  # n_events (out)
        ctypes.c_uint32,  # window_length1
        ctypes.c_uint32,  # window_length2
        ctypes.c_float,  # threshold1
        ctypes.c_float,  # threshold2
        ctypes.c_float,  # peak_height
        ctypes.c_uint32,  # min_seg_len
        ctypes.c_uint32,  # max_seg_len
        ctypes.POINTER(ctypes.c_double),  # mean_sum (in/out)
        ctypes.POINTER(ctypes.c_double),  # std_dev_sum (in/out)
        ctypes.POINTER(ctypes.c_uint32),  # n_events_sum (in/out)
    ]
    _RH2_LIB.rh2_free.argtypes = [ctypes.c_void_p]
    _RH2_LIB.rh2_free.restype = None

    return _RH2_LIB


# ── RawHash v1.0 defaults (from extern/RawHash/src/roptions.c) ──────────────

RH1_DEFAULTS = {
    "window_length1": 3,
    "window_length2": 6,
    "threshold1": 4.30265,
    "threshold2": 2.57058,
    "peak_height": 1.0,
}

# ── RawHash2 defaults (from extern/RawHash2/src/roptions.c) ─────────────────

RH2_DEFAULTS = {
    "window_length1": 3,
    "window_length2": 9,
    "threshold1": 4.0,
    "threshold2": 3.5,
    "peak_height": 0.4,
    "min_segment_length": 0,
    "max_segment_length": 500,
}


# ── Public API ───────────────────────────────────────────────────────────────


def detect_events_rawhash(
    signal_pA: np.ndarray,
    window_length1: int = RH1_DEFAULTS["window_length1"],
    window_length2: int = RH1_DEFAULTS["window_length2"],
    threshold1: float = RH1_DEFAULTS["threshold1"],
    threshold2: float = RH1_DEFAULTS["threshold2"],
    peak_height: float = RH1_DEFAULTS["peak_height"],
) -> np.ndarray:
    """
    Run RawHash v1.0 event detection on a pA signal.

    The signal should already be calibrated to pA and filtered to 30-200 range
    (this is done by signal_reader.py).

    Returns a 1D float32 array of event values (z-score normalized).
    """
    lib = _get_rh1_lib()
    sig = np.ascontiguousarray(signal_pA, dtype=np.float32)
    n_events = ctypes.c_uint32(0)

    sig_ptr = sig.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    result_ptr = lib.rh1_detect_events(
        ctypes.c_uint32(len(sig)),
        sig_ptr,
        ctypes.byref(n_events),
        ctypes.c_uint32(window_length1),
        ctypes.c_uint32(window_length2),
        ctypes.c_float(threshold1),
        ctypes.c_float(threshold2),
        ctypes.c_float(peak_height),
    )

    n = n_events.value
    if not result_ptr or n == 0:
        if result_ptr:
            lib.rh1_free(result_ptr)
        return np.array([], dtype=np.float32)

    events = np.ctypeslib.as_array(result_ptr, shape=(n,)).copy()
    lib.rh1_free(result_ptr)
    return events


def detect_events_rawhash2(
    signal_pA: np.ndarray,
    window_length1: int = RH2_DEFAULTS["window_length1"],
    window_length2: int = RH2_DEFAULTS["window_length2"],
    threshold1: float = RH2_DEFAULTS["threshold1"],
    threshold2: float = RH2_DEFAULTS["threshold2"],
    peak_height: float = RH2_DEFAULTS["peak_height"],
    min_segment_length: int = RH2_DEFAULTS["min_segment_length"],
    max_segment_length: int = RH2_DEFAULTS["max_segment_length"],
) -> np.ndarray:
    """
    Run RawHash2 event detection on a pA signal.

    The signal should already be calibrated to pA and filtered to 30-200 range.

    Returns a 1D float32 array of event values.
    """
    lib = _get_rh2_lib()
    sig = np.ascontiguousarray(signal_pA, dtype=np.float32)
    n_events = ctypes.c_uint32(0)

    c_mean_sum = ctypes.c_double(0.0)
    c_std_dev_sum = ctypes.c_double(0.0)
    c_n_events_sum = ctypes.c_uint32(0)

    sig_ptr = sig.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    result_ptr = lib.rh2_detect_events(
        ctypes.c_uint32(len(sig)),
        sig_ptr,
        ctypes.byref(n_events),
        ctypes.c_uint32(window_length1),
        ctypes.c_uint32(window_length2),
        ctypes.c_float(threshold1),
        ctypes.c_float(threshold2),
        ctypes.c_float(peak_height),
        ctypes.c_uint32(min_segment_length),
        ctypes.c_uint32(max_segment_length),
        ctypes.byref(c_mean_sum),
        ctypes.byref(c_std_dev_sum),
        ctypes.byref(c_n_events_sum),
    )

    n = n_events.value
    if not result_ptr or n == 0:
        if result_ptr:
            lib.rh2_free(result_ptr)
        return np.array([], dtype=np.float32)

    events = np.ctypeslib.as_array(result_ptr, shape=(n,)).copy()
    lib.rh2_free(result_ptr)
    return events
