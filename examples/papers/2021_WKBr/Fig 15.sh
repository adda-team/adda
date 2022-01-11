#!/bin/bash

# This is an example. 
# We do not use checkpoints. It is assumed that computer / cluster resources for calculations are provided for an unlimited time. 
# Otherwise, you need to use checkpoints (read the manual).
# WKBR_I_PATH, WKBR_II_PATH - the paths to the files of internal electric fields, calculated using the WKBr method (the case of a sphere).
# EXACT_PATH - the path to the file of exact internal electric fields, calculated using Bhfield (see necessary links in our paper).

# Fig.15 [a, b, c]
adda_mpi.exe -shape sphere -size 100 -grid 180 -m [1.01, 1.05, 1.1] 0 -no_vol_cor -eps 6 -init_field read $EXACT_PATH
adda_mpi.exe -shape sphere -size 100 -grid 180 -m [1.01, 1.05, 1.1] 0 -no_vol_cor -eps 6 -init_field inc
adda_mpi.exe -shape sphere -size 100 -grid 180 -m [1.01, 1.05, 1.1] 0 -no_vol_cor -eps 6 -init_field zero
adda_mpi.exe -shape sphere -size 100 -grid 180 -m [1.01, 1.05, 1.1] 0 -no_vol_cor -eps 6 -init_field wkb
adda_mpi.exe -shape sphere -size 100 -grid 180 -m [1.01, 1.05, 1.1] 0 -no_vol_cor -eps 6 -init_field read $WKBR_I_PATH
adda_mpi.exe -shape sphere -size 100 -grid 180 -m [1.01, 1.05, 1.1] 0 -no_vol_cor -eps 6 -init_field read $WKBR_II_PATH