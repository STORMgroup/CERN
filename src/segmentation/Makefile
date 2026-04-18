# Thin Makefile wrapper around the CMake build.
#
# Options (pass-through to CMake):
#   CMAKE_OPTS="-DENABLE_HDF5=ON -DENABLE_SLOW5=ON"
# Examples:
#   make                                                 # POD5 only (default)
#   make CMAKE_OPTS="-DENABLE_HDF5=ON -DENABLE_SLOW5=ON" # all formats
#   make clean

.PHONY: all cmake clean help

ROOT_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
JOBS     ?= $(shell nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
CMAKE_OPTS ?=

all: cmake

cmake:
	@mkdir -p build
	cd build && cmake $(CMAKE_OPTS) $(ROOT_DIR)
	cmake --build build -j$(JOBS)
	@mkdir -p bin && cp build/src/generate_events bin/
	@echo "Build complete: bin/generate_events"

clean:
	rm -rf bin build

help:
	@echo "Segmentation build system (CMake wrapper)"
	@echo ""
	@echo "Targets:"
	@echo "  make                                          Build with POD5 enabled (default)"
	@echo '  make CMAKE_OPTS="-DENABLE_HDF5=ON"            Enable FAST5 support'
	@echo '  make CMAKE_OPTS="-DENABLE_SLOW5=ON"           Enable SLOW5/BLOW5 support'
	@echo '  make CMAKE_OPTS="-DENABLE_POD5=OFF ..."       Disable POD5'
	@echo "  make clean                                    Remove build/ and bin/"
