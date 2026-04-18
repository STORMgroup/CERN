# Dependencies.cmake
# Orchestrates building/fetching all external dependencies based on feature toggles.
#
# Feature toggles (set in root CMakeLists.txt):
#   ENABLE_HDF5   - build with HDF5/FAST5 support
#   ENABLE_SLOW5  - build with SLOW5/BLOW5 support
#   ENABLE_POD5   - build with POD5 support

# --- Submodule check ---
function(check_submodule path name)
    if(NOT EXISTS "${path}/.git" AND NOT EXISTS "${path}/CMakeLists.txt" AND NOT EXISTS "${path}/Makefile")
        message(FATAL_ERROR
            "Git submodule '${name}' not found at ${path}\n"
            "Please initialize submodules:\n"
            "  git submodule update --init --recursive\n")
    endif()
endfunction()

# --- zstd (needed by POD5 and slow5lib for BLOW5 compression) ---
# zstd MUST be built before slow5lib so that slow5lib can find zstd headers.
if(ENABLE_POD5 OR ENABLE_SLOW5)
    set(ZSTD_SOURCE_DIR "${CMAKE_SOURCE_DIR}/extern/zstd")
    check_submodule("${ZSTD_SOURCE_DIR}" "zstd")

    if(EXISTS "${ZSTD_SOURCE_DIR}/build/cmake/CMakeLists.txt")
        # Older zstd versions declare cmake_minimum_required < 3.5 which CMake 4.x
        # rejects; CMAKE_POLICY_VERSION_MINIMUM lets them configure anyway.
        set(CMAKE_POLICY_VERSION_MINIMUM 3.5 CACHE STRING "" FORCE)
        set(ZSTD_BUILD_SHARED   OFF CACHE BOOL "" FORCE)
        set(ZSTD_BUILD_PROGRAMS OFF CACHE BOOL "" FORCE)
        set(ZSTD_BUILD_TESTS    OFF CACHE BOOL "" FORCE)
        add_subdirectory("${ZSTD_SOURCE_DIR}/build/cmake" "${CMAKE_BINARY_DIR}/extern/zstd" EXCLUDE_FROM_ALL)
        set(ZSTD_FOUND TRUE)
        set(ZSTD_LIB_DIR "${CMAKE_BINARY_DIR}/extern/zstd/lib")
        message(STATUS "zstd will be built from submodule (CMake): ${ZSTD_SOURCE_DIR}")
    else()
        include(ExternalProject)
        include(ProcessorCount)
        ProcessorCount(NPROC)
        if(NPROC EQUAL 0)
            set(NPROC 1)
        endif()
        ExternalProject_Add(zstd_external
            SOURCE_DIR      "${ZSTD_SOURCE_DIR}"
            CONFIGURE_COMMAND ""
            BUILD_COMMAND     make -j${NPROC}
            INSTALL_COMMAND   ""
            BUILD_IN_SOURCE   TRUE
            LOG_BUILD         TRUE
        )
        set(ZSTD_FOUND TRUE)
        set(ZSTD_LIB_DIR "${ZSTD_SOURCE_DIR}/lib")
        message(STATUS "zstd will be built from submodule (Make fallback): ${ZSTD_SOURCE_DIR}")
    endif()
endif()

# --- HDF5 ---
if(ENABLE_HDF5)
    if(NOT USE_SYSTEM_HDF5)
        check_submodule("${CMAKE_SOURCE_DIR}/extern/hdf5" "hdf5")
    endif()
    include(cmake/FetchHDF5.cmake)
endif()

# --- slow5lib ---
if(ENABLE_SLOW5)
    check_submodule("${CMAKE_SOURCE_DIR}/extern/slow5lib" "slow5lib")
    include(cmake/FetchSlow5.cmake)
endif()

# --- POD5 (pre-built download, no submodule) ---
if(ENABLE_POD5)
    include(cmake/FetchPOD5.cmake)
endif()
