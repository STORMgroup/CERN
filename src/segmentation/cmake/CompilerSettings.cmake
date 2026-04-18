# CompilerSettings.cmake
# Optimization flags, sanitizers, and warnings.

include(CheckCXXCompilerFlag)

message(STATUS "System: ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")

# --- Optimization flags ---
check_cxx_compiler_flag("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(ARCH_FLAG "-march=native")
else()
    set(ARCH_FLAG "")
    message(STATUS "Compiler does not support -march=native, using default architecture flags")
endif()

set(CMAKE_C_FLAGS_RELEASE   "-O3 ${ARCH_FLAG}" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 ${ARCH_FLAG}" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_DEBUG     "-O2 -g ${ARCH_FLAG}" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_DEBUG   "-O2 -g ${ARCH_FLAG}" CACHE STRING "" FORCE)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    message(STATUS "Build type not specified, defaulting to Release")
endif()

# --- Sanitizer support ---
option(ENABLE_ASAN "Enable AddressSanitizer" OFF)
option(ENABLE_TSAN "Enable ThreadSanitizer" OFF)

if(ENABLE_ASAN)
    add_compile_options(-fsanitize=address)
    add_link_options(-fsanitize=address)
    message(STATUS "AddressSanitizer enabled")
endif()
if(ENABLE_TSAN)
    add_compile_options(-fsanitize=thread)
    add_link_options(-fsanitize=thread)
    message(STATUS "ThreadSanitizer enabled")
endif()

add_compile_options(-Wall)
