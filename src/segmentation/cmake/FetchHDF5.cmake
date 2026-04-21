# FetchHDF5.cmake
# Build HDF5 from the extern/hdf5 git submodule, or use a system-installed HDF5.
#
# Sets:
#   HDF5_INCLUDE_DIRS  - include path for HDF5 headers
#   HDF5_LIBRARIES     - static library path for HDF5
#   HDF5_FOUND         - TRUE if HDF5 is available

include(ExternalProject)

option(USE_SYSTEM_HDF5 "Use system-installed HDF5 instead of building from submodule" OFF)

if(USE_SYSTEM_HDF5)
    find_package(HDF5 REQUIRED COMPONENTS C)
    if(HDF5_FOUND)
        message(STATUS "Using system HDF5: ${HDF5_INCLUDE_DIRS}")
    endif()
else()
    set(HDF5_SOURCE_DIR "${CMAKE_SOURCE_DIR}/extern/hdf5")
    set(HDF5_INSTALL_DIR "${CMAKE_BINARY_DIR}/extern/hdf5-install")
    set(HDF5_INCLUDE_DIRS "${HDF5_INSTALL_DIR}/include")
    set(HDF5_LIBRARIES "${HDF5_INSTALL_DIR}/lib/libhdf5.a")

    include(ProcessorCount)
    ProcessorCount(NPROC)
    if(NPROC EQUAL 0)
        set(NPROC 1)
    endif()

    ExternalProject_Add(hdf5_external
        SOURCE_DIR      "${HDF5_SOURCE_DIR}"
        INSTALL_DIR     "${HDF5_INSTALL_DIR}"
        CONFIGURE_COMMAND
            <SOURCE_DIR>/configure
                --enable-threadsafe
                --disable-hl
                --prefix=<INSTALL_DIR>
        BUILD_COMMAND   make -j${NPROC}
        INSTALL_COMMAND make install
        BUILD_BYPRODUCTS
            "${HDF5_INSTALL_DIR}/lib/libhdf5.a"
            "${HDF5_INSTALL_DIR}/lib64/libhdf5.a"
        LOG_CONFIGURE   TRUE
        LOG_BUILD       TRUE
        LOG_INSTALL     TRUE
    )

    set(HDF5_FOUND TRUE)
    message(STATUS "HDF5 will be built from submodule: ${HDF5_SOURCE_DIR}")
endif()
