﻿# CMakeList.txt : CMake project for SpectralEntropy_C, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

project ("SpectralEntropy")

# Add source to this project's executable.
ADD_LIBRARY(SpectralEntropy SHARED "SpectralEntropy.c" "SpectralEntropy.h" "CleanSpectrum.c" "CleanSpectrum.h")
add_executable (SpectralEntropyTest "Main.c")
TARGET_LINK_LIBRARIES(SpectralEntropyTest SpectralEntropy m)


if (CMAKE_VERSION VERSION_GREATER 3.12)
    set_property(TARGET SpectralEntropyTest PROPERTY C_STANDARD 99)
endif()
