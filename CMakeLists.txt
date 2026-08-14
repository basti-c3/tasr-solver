cmake_minimum_required(VERSION 3.20)

project(Solver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(ENABLE_HIGHS "Compile the HiGHS backend" ON)
option(ENABLE_GUROBI "Compile the Gurobi C-API backend" ON)

if(NOT ENABLE_HIGHS AND NOT ENABLE_GUROBI)
    message(FATAL_ERROR "Enable at least one solver backend")
endif()

# -------------------------------------------------
# Make packages installed inside the active MSYS2
# environment discoverable by native CMake.
# Example compiler:
#   C:/msys64/ucrt64/bin/c++.exe
# -> prefix:
#   C:/msys64/ucrt64
# -------------------------------------------------
if(MINGW)
    get_filename_component(_compiler_bin "${CMAKE_CXX_COMPILER}" DIRECTORY)
    get_filename_component(_mingw_prefix "${_compiler_bin}" DIRECTORY)
    list(PREPEND CMAKE_PREFIX_PATH "${_mingw_prefix}")
endif()

# -------------------------------------------------
# Debug configuration
# -------------------------------------------------
set(DEBUG_INSTANCE "setA-01")
set(DEBUG_ITERATIONS 1)
set(DEBUG_TIME_LIMIT 600.0)
set(DEBUG_SOLVER "highs")

set(DEBUG_OUTPUT_FILE "instances/${DEBUG_INSTANCE}-srpaths.json")
set(DEBUG_MODEL_PREFIX "debug-${DEBUG_INSTANCE}")

# -------------------------------------------------
# Solver target
# -------------------------------------------------
add_executable(
    solver
    src/main.cpp
    src/instance_data.cpp
    src/split_coefficients.cpp
    src/tasr_solver.cpp
    src/solver_backend.cpp
)

target_compile_definitions(
    solver
    PRIVATE
    NOMINMAX
    _USE_MATH_DEFINES
    _SH_DENYNO=0

    $<$<CONFIG:Debug>:WRITE_LP_MODELS>
    $<$<CONFIG:Debug>:DEBUG_DEFAULT_INSTANCE="${DEBUG_INSTANCE}">
    $<$<CONFIG:Debug>:DEBUG_DEFAULT_OUTPUT_FILE="${DEBUG_OUTPUT_FILE}">
    $<$<CONFIG:Debug>:DEBUG_DEFAULT_ITERATIONS=${DEBUG_ITERATIONS}>
    $<$<CONFIG:Debug>:DEBUG_DEFAULT_TIME_LIMIT=${DEBUG_TIME_LIMIT}>
    $<$<CONFIG:Debug>:DEBUG_DEFAULT_SOLVER="${DEBUG_SOLVER}">
    $<$<CONFIG:Debug>:DEBUG_DEFAULT_MODEL_PREFIX="${DEBUG_MODEL_PREFIX}">
)

target_include_directories(
    solver
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/external/networktools/networktools
)

# -------------------------------------------------
# HiGHS backend
# -------------------------------------------------
if(ENABLE_HIGHS)
    find_package(highs CONFIG REQUIRED)
    target_link_libraries(solver PRIVATE highs::highs)
    target_compile_definitions(solver PRIVATE ROADEF_HAVE_HIGHS=1)
    message(STATUS "HiGHS backend enabled")
endif()

# -------------------------------------------------
# Gurobi backend - C API only
#
# This deliberately does NOT use gurobi_c++.h or the Gurobi C++ library.
# The C-API DLL boundary is what we tested successfully from MSYS2 UCRT64.
# -------------------------------------------------
if(ENABLE_GUROBI)
    set(GUROBI_ROOT "$ENV{GUROBI_HOME}" CACHE PATH "Gurobi installation root")

    if(NOT GUROBI_ROOT)
        message(FATAL_ERROR
            "ENABLE_GUROBI=ON but GUROBI_HOME is not set. "
            "Set GUROBI_HOME in your MSYS2 ~/.bashrc or configure with "
            "-DENABLE_GUROBI=OFF.")
    endif()

    # If GUROBI_HOME is stored in MSYS form such as /c/gurobi1302/win64,
    # convert it to a Windows/CMake path for native CMake.
    if(MINGW AND GUROBI_ROOT MATCHES "^/")
        find_program(CYGPATH_EXECUTABLE cygpath)
        if(CYGPATH_EXECUTABLE)
            execute_process(
                COMMAND "${CYGPATH_EXECUTABLE}" -m "${GUROBI_ROOT}"
                OUTPUT_VARIABLE _gurobi_windows_path
                OUTPUT_STRIP_TRAILING_WHITESPACE
                RESULT_VARIABLE _cygpath_result
            )
            if(_cygpath_result EQUAL 0 AND _gurobi_windows_path)
                set(GUROBI_ROOT "${_gurobi_windows_path}")
            endif()
        endif()
    endif()

    message(STATUS "GUROBI_ROOT = ${GUROBI_ROOT}")

    find_path(
        GUROBI_INCLUDE_DIR
        NAMES gurobi_c.h
        HINTS "${GUROBI_ROOT}"
        PATH_SUFFIXES include
        REQUIRED
    )

    # Gurobi's versioned core library. Add future versions here if needed.
    find_library(
        GUROBI_C_LIBRARY
        NAMES
            gurobi140
            gurobi130
            gurobi120
            gurobi110
            gurobi100
            gurobi95
            gurobi
        HINTS "${GUROBI_ROOT}"
        PATH_SUFFIXES lib bin
        REQUIRED
    )

    target_include_directories(solver PRIVATE "${GUROBI_INCLUDE_DIR}")
    target_link_libraries(solver PRIVATE "${GUROBI_C_LIBRARY}")
    target_compile_definitions(solver PRIVATE ROADEF_HAVE_GUROBI=1)

    message(STATUS "Gurobi C-API backend enabled")
    message(STATUS "Gurobi include: ${GUROBI_INCLUDE_DIR}")
    message(STATUS "Gurobi library: ${GUROBI_C_LIBRARY}")
endif()
