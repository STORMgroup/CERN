# FetchPOD5.cmake
# Download pre-built POD5 libraries for the current platform.
#
# Configurable:
#   POD5_VERSION     - version to download (default: 0.3.36)
#
# Sets:
#   POD5_INCLUDE_DIRS  - include path for POD5 headers
#   POD5_LIBRARIES     - list of static libraries to link
#   POD5_FOUND         - TRUE if POD5 is available

include(ExternalProject)

set(POD5_REPO "https://github.com/nanoporetech/pod5-file-format")
set(POD5_DOWNLOAD_DIR "${CMAKE_BINARY_DIR}/extern/pod5-${POD5_VERSION}")

if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm64|ARM64")
        set(POD5_URL "${POD5_REPO}/releases/download/${POD5_VERSION}/lib_pod5-${POD5_VERSION}-linux-arm64.tar.gz")
    else()
        set(POD5_URL "${POD5_REPO}/releases/download/${POD5_VERSION}/lib_pod5-${POD5_VERSION}-linux-x64.tar.gz")
    endif()
    set(POD5_LIB_SUBDIR "lib64")
    set(POD5_EXTRA_LIBS "${POD5_DOWNLOAD_DIR}/lib64/libjemalloc_pic.a")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|ARM64|aarch64")
        set(POD5_URL "${POD5_REPO}/releases/download/${POD5_VERSION}/lib_pod5-${POD5_VERSION}-osx-14.0-arm64.tar.gz")
    else()
        message(FATAL_ERROR "POD5 ${POD5_VERSION} does not provide pre-built binaries for macOS x86_64. "
                            "Use -DENABLE_POD5=OFF to build without POD5 support.")
    endif()
    set(POD5_LIB_SUBDIR "lib")
    set(POD5_EXTRA_LIBS "")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    set(POD5_URL "${POD5_REPO}/releases/download/${POD5_VERSION}/lib_pod5-${POD5_VERSION}-win-x64.tar.gz")
    set(POD5_LIB_SUBDIR "lib")
    set(POD5_EXTRA_LIBS "")
else()
    message(FATAL_ERROR "Unsupported platform for POD5: ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}")
endif()

set(POD5_INCLUDE_DIRS "${POD5_DOWNLOAD_DIR}/include")
set(POD5_LIB_DIR "${POD5_DOWNLOAD_DIR}/${POD5_LIB_SUBDIR}")
set(POD5_LIBRARIES
    "${POD5_LIB_DIR}/libpod5_format.a"
    "${POD5_LIB_DIR}/libarrow.a"
    ${POD5_EXTRA_LIBS}
)

message(STATUS "POD5 v${POD5_VERSION} will be downloaded from: ${POD5_URL}")

ExternalProject_Add(pod5_external
    URL             "${POD5_URL}"
    SOURCE_DIR      "${POD5_DOWNLOAD_DIR}"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND     ""
    INSTALL_COMMAND   ""
    BUILD_BYPRODUCTS  ${POD5_LIBRARIES}
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    LOG_DOWNLOAD      TRUE
)

set(POD5_FOUND TRUE)
