# src/CMakeLists.txt — generate_events binary

set(SEG_SOURCES
    rawhash_events.c
    rawhash2_events.c
    signal_reader.cpp
    main.cpp
)

# Per-format readers are added to the source list only when their feature is
# enabled; otherwise signal_reader.cpp provides a stub that throws at runtime.
if(ENABLE_POD5)
    list(APPEND SEG_SOURCES pod5_reader.cpp)
endif()
if(ENABLE_HDF5)
    list(APPEND SEG_SOURCES fast5_reader.cpp)
endif()
if(ENABLE_SLOW5)
    list(APPEND SEG_SOURCES slow5_reader.cpp)
endif()

add_executable(generate_events ${SEG_SOURCES})

# POD5's pre-built archives need the C++ linker even if main.cpp wouldn't
# otherwise force it. Same trick RawHash2 uses: compile the C sources as C++
# when at least one C++ backend is enabled. Without any C++ backend, the C
# files stay C and link cleanly against libm.
if(ENABLE_POD5 OR ENABLE_HDF5)
    set_source_files_properties(rawhash_events.c rawhash2_events.c
                                PROPERTIES LANGUAGE CXX)
endif()

target_include_directories(generate_events PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# --- POD5 ------------------------------------------------------------------
if(ENABLE_POD5)
    target_compile_definitions(generate_events PRIVATE SEG_HAVE_POD5=1)
    target_include_directories(generate_events PRIVATE ${POD5_INCLUDE_DIRS})
    target_link_libraries(generate_events PRIVATE ${POD5_LIBRARIES})
    add_dependencies(generate_events pod5_external)
endif()

# --- HDF5 ------------------------------------------------------------------
if(ENABLE_HDF5)
    target_compile_definitions(generate_events PRIVATE SEG_HAVE_HDF5=1)
    target_include_directories(generate_events PRIVATE ${HDF5_INCLUDE_DIRS})
    target_link_libraries(generate_events PRIVATE ${HDF5_LIBRARIES})
    if(NOT USE_SYSTEM_HDF5 AND TARGET hdf5_external)
        add_dependencies(generate_events hdf5_external)
    endif()
endif()

# --- SLOW5 -----------------------------------------------------------------
if(ENABLE_SLOW5)
    target_compile_definitions(generate_events PRIVATE SEG_HAVE_SLOW5=1)
    target_include_directories(generate_events PRIVATE ${SLOW5_INCLUDE_DIRS})
    target_link_libraries(generate_events PRIVATE ${SLOW5_LIBRARIES})
    add_dependencies(generate_events slow5_external)
endif()

# --- zstd (required by POD5 and SLOW5 for BLOW5 decompression) -------------
if(ENABLE_POD5 OR ENABLE_SLOW5)
    if(TARGET libzstd_static)
        target_link_libraries(generate_events PRIVATE libzstd_static)
    elseif(ZSTD_LIB_DIR)
        target_link_libraries(generate_events PRIVATE "${ZSTD_LIB_DIR}/libzstd.a")
        if(TARGET zstd_external)
            add_dependencies(generate_events zstd_external)
        endif()
    endif()
endif()

# --- System libraries ------------------------------------------------------
find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)

target_link_libraries(generate_events PRIVATE
    Threads::Threads
    ZLIB::ZLIB
    ${CMAKE_DL_LIBS}
    m
)

# POD5 static archives reference C++ stdlib symbols that must resolve after
# the archives on Linux (GCC < 9: -lstdc++fs for std::filesystem; always:
# trailing -lstdc++ for correct archive → runtime link order).
if(ENABLE_POD5 OR ENABLE_HDF5)
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        target_link_libraries(generate_events PRIVATE stdc++fs stdc++)
    endif()
endif()

install(TARGETS generate_events RUNTIME DESTINATION bin)
