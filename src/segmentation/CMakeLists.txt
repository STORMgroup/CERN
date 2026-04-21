cmake_minimum_required(VERSION 3.16)
project(segmentation
    VERSION 1.0
    LANGUAGES C CXX
    DESCRIPTION "Nanopore signal event segmentation (RawHash / RawHash2)"
)

# --- Language standards ---
# C++17 is needed for std::filesystem (used by the directory-iteration reader).
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD 99)

# --- Feature toggles ---
# Control which nanopore input formats the binary will support.
# Usage: cmake -DENABLE_HDF5=ON -DENABLE_SLOW5=ON ..
option(ENABLE_HDF5  "Enable HDF5/FAST5 support"      OFF)
option(ENABLE_SLOW5 "Enable SLOW5/BLOW5 support"     OFF)
option(ENABLE_POD5  "Enable POD5 support"            ON)

# --- Dependency version parameters ---
# Override at configure time, e.g. `cmake -DPOD5_VERSION=0.3.11 ..`
set(POD5_VERSION "0.3.36" CACHE STRING "POD5 library version to download")

# --- Compiler and platform settings ---
include(cmake/CompilerSettings.cmake)

# --- External dependencies ---
include(cmake/Dependencies.cmake)

# --- Target ---
add_subdirectory(src)
