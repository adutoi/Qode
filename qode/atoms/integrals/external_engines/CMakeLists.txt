cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

# provide some sort of BLAS or LAPACK via e.g. IntelMKL
#
# source /opt/intel/oneapi/mkl/latest/env/vars.sh
# export MKL_DIR=${MKLROOT}
#
# then build into a build directory and the following setup
# will take care of the rest
 
project(pylible_ints
    LANGUAGES CXX
    DESCRIPTION "Python bindings for the lible AO integral library")
 
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
 
# ---------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------
 
# FetchContent is a standard CMake module; it must be included before use.
include(FetchContent)
 
# pybind11 – find via cmake find_package (pybind11 must be installed, e.g.
# via `pip install pybind11` or a system package).
find_package(pybind11 REQUIRED)
 
# OpenMP (lible uses it internally; linking against it avoids undefined refs
# at module load time).
find_package(OpenMP REQUIRED)
 
# lible – fetched from GitHub at configure time.
# Pin to a specific commit for reproducibility; update the GIT_TAG as needed.
FetchContent_Declare(lible
    GIT_REPOSITORY https://github.com/MihkuU/Lible
    GIT_TAG        d277c98
)
FetchContent_MakeAvailable(lible)
 
# ---------------------------------------------------------------------------
# Module target
# ---------------------------------------------------------------------------
 
# pybind11_add_module must come BEFORE any target_link_libraries calls that
# reference the target it creates (pylible_ints).
pybind11_add_module(pylible_ints
    ${CMAKE_CURRENT_SOURCE_DIR}/pylible_ints.cpp
)
 
target_link_libraries(pylible_ints PRIVATE
    lible::lible
    OpenMP::OpenMP_CXX
)
 
# Ensure C++20 is used consistently with lible.
target_compile_features(pylible_ints PRIVATE cxx_std_20)

# put targets into same directory
set_target_properties(pylible_ints PROPERTIES
    BUILD_RPATH "\$ORIGIN"
    INSTALL_RPATH "\$ORIGIN"
)
 
# ---------------------------------------------------------------------------
# Install
#
# The built .so must be importable by Python as `import pylible_ints`.
# By default we install it into the external_engines source directory so that
# the existing `import pylible_ints as pli` in pybind_hummr_ints.py works
# without any PYTHONPATH changes.
#
# Override at configure time with:
#   cmake -DPYLIBLE_INSTALL_DIR=/your/site-packages ...
# ---------------------------------------------------------------------------
 
if(NOT DEFINED PYLIBLE_INSTALL_DIR)
    set(PYLIBLE_INSTALL_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
endif()
 
install(TARGETS pylible_ints
    LIBRARY DESTINATION "${PYLIBLE_INSTALL_DIR}"
)

install(TARGETS lible
    LIBRARY DESTINATION "${PYLIBLE_INSTALL_DIR}"
)
