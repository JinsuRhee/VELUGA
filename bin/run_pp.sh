#!/bin/bash

# -----------------------------
#  keyword_input
# -----------------------------
read -p "Enter headername: " fname
read -p "Enter horg (h or g): " horg
read -p "Num_thread?: " nth

# check 
if [[ "$horg" != "h" && "$horg" != "g" ]]; then
    echo "Error: horg must be 'h' or 'g'"
    exit 1
fi

# -----------------------------
# HDF5 library
# -----------------------------
if [[ -z "$CONDA_PREFIX" ]]; then
    echo "Error: CONDA_PREFIX is not set. Please activate your conda environment."
    exit 1
fi

HDF5_LIB="$CONDA_PREFIX/lib/libhdf5.so"

# -----------------------------
# IDL
# -----------------------------
LD_PRELOAD=$HDF5_LIB idl -e "run_pp, '$fname', horg='$horg', num_thread=$nth"

