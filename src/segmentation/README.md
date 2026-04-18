# Segmentation — Event Generation

Generate precomputed event values from nanopore signal files using the
RawHash (ScrappieR9) or RawHash2 (ScrappieR10) segmentation algorithms.

This is the C++/CMake port of the original Python wrapper; the legacy Python
implementation is preserved in [`python-legacy/`](python-legacy/) for
reference.

## Build

```bash
# From the repo root
git submodule update --init --recursive

cd src/segmentation
make                                                   # POD5 only (default)
make CMAKE_OPTS="-DENABLE_HDF5=ON -DENABLE_SLOW5=ON"   # all formats
```

The binary lands at `src/segmentation/bin/generate_events`.

### Feature toggles

The build is modular — each input format is opt-in at CMake configure time.
A binary built without a given format will report a runtime error when given
that extension.

| Flag             | Default | Enables    |
|------------------|---------|------------|
| `-DENABLE_POD5`  | `ON`    | `.pod5`    |
| `-DENABLE_HDF5`  | `OFF`   | `.fast5`   |
| `-DENABLE_SLOW5` | `OFF`   | `.slow5`, `.blow5` |

POD5 is downloaded as a pre-built archive (default `POD5_VERSION=0.3.36`);
HDF5 and SLOW5 are built from the submodules under `extern/`.

### Requirements

- CMake ≥ 3.16, a C++17 compiler (GCC 9+ / Clang 10+), `make`, `git`.
- For the HDF5 submodule build: autotools (`autoconf`, `automake`, `libtool`).
- For SLOW5: its submodule pulls in `zstd` (already a submodule here).

If you prefer a system-installed HDF5 (e.g. `libhdf5-dev`), configure with
`-DUSE_SYSTEM_HDF5=ON` to skip the submodule build.

## Usage

```bash
# RawHash2 (ScrappieR10) with R10.4.1 parameters
bin/generate_events -m rawhash2 --r10 -i /path/to/signals.pod5 -o events.tsv

# RawHash2 with base defaults (without --r10 preset)
bin/generate_events -m rawhash2 -i /path/to/signals.pod5 -o events.tsv

# RawHash v1.0 (ScrappieR9) on a directory of signal files
bin/generate_events -m rawhash -i /path/to/signal_dir/ -o events.tsv
```

When given a directory, all signal files in that directory are processed.
Run `bin/generate_events --help` for the full list of parameter overrides.

### Supported input formats

| Format | Extensions         | Library       |
|--------|--------------------|---------------|
| FAST5  | `.fast5`           | libhdf5       |
| POD5   | `.pod5`            | pod5 (C API)  |
| SLOW5  | `.slow5`, `.blow5` | slow5lib      |

### Output format

Tab-separated values (TSV), one read per line, floats formatted as `%.6f`:

```
read_id1\tevent1\tevent2\t…\teventn
read_id2\tevent1\tevent2\t…\teventn
```

Byte-for-byte compatible with the legacy Python output and with the
`--events-file` flag in [RawHash2](https://github.com/STORMgroup/RawHash2).

## Algorithm details

### RawHash v1.0 (`--method rawhash`)

1. **Input**: Raw signal calibrated to pA, filtered to 30–200 pA.
2. **T-test segmentation**: Two sliding windows detect change-points on the
   pA signal directly (no pre-normalization).
3. **Event computation**: Simple mean of signal samples between consecutive
   peaks (using prefix sums).
4. **Post-normalization**: Z-score normalization of the generated events.

Default parameters (from `extern/RawHash/src/roptions.c`):
- `window_length1=3`, `threshold1=4.30265`
- `window_length2=6`, `threshold2=2.57058`
- `peak_height=1.0`

### RawHash2 (`--method rawhash2`)

1. **Input**: Raw signal calibrated to pA, filtered to 30–200 pA.
2. **Z-score normalization**: Normalize pA signal, then discard samples
   outside ±3σ.
3. **T-test segmentation**: Two sliding windows detect change-points on the
   normalized signal.
4. **Event computation**: IQR-filtered robust mean of normalized signal
   between peaks, with min/max segment length filtering.

Default parameters (from `extern/RawHash2/src/roptions.c`):
- `window_length1=3`, `threshold1=4.0`
- `window_length2=9`, `threshold2=3.5`
- `peak_height=0.4`
- `min_segment_length=0`, `max_segment_length=500`

**R10.4.1 preset** (`--r10`):
- `window_length1=3`, `threshold1=6.5`
- `window_length2=6`, `threshold2=4.0`
- `peak_height=0.2`
- `min_segment_length=0`, `max_segment_length=500`

Use `--r10` when generating events for R10.4.1 data to match the parameters
used by `rawhash2 --r10`.

## Structure

```
src/segmentation/
├── CMakeLists.txt              # Top-level CMake
├── Makefile                    # Thin `make cmake` wrapper
├── cmake/
│   ├── CompilerSettings.cmake  # -O3 -march=native, sanitizers, warnings
│   ├── Dependencies.cmake      # Feature-toggled dependency orchestration
│   ├── FetchPOD5.cmake         # Pre-built POD5 download
│   ├── FetchSlow5.cmake        # slow5lib submodule build
│   └── FetchHDF5.cmake         # HDF5 submodule (or system) build
├── extern/
│   ├── slow5lib/               # submodule — hasindu2008/slow5lib
│   ├── hdf5/                   # submodule — HDFGroup/hdf5
│   └── zstd/                   # submodule — facebook/zstd
├── src/
│   ├── main.cpp                # CLI entry point
│   ├── rawhash_events.{c,h}    # RawHash v1.0 event detection
│   ├── rawhash2_events.{c,h}   # RawHash2 event detection
│   ├── signal_reader.{hpp,cpp} # Extension dispatcher + 30–200 pA filter
│   ├── fast5_reader.cpp        # HDF5 C API — compiled iff ENABLE_HDF5
│   ├── pod5_reader.cpp         # POD5 C API — compiled iff ENABLE_POD5
│   └── slow5_reader.cpp        # slow5lib      — compiled iff ENABLE_SLOW5
├── python-legacy/              # Preserved Python+ctypes implementation
└── bin/generate_events         # Build artifact
```
