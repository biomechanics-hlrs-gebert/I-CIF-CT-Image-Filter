#!/bin/bash
# ----------------------------------------------------------------------------------------
# Set Environment for Image Processing
#
# Author:          Johannes Gebert - HLRS - NUM «gebert@hlrs.de»
# Created:         25.04.2021
# Last edit:       25.04.2021
#
# Usage: source «this script»
# ----------------------------------------------------------------------------------------
#
# Set Architecture
# hawk    - HPE Apollo
# vulcan  - NEC Cluster
# julius  - A Whiskey Lake Notebook, 4 cores, 16Gb memory, APU
export IP_ARCH="julius"
export IP_VRSN="V0.2.0"
# ----------------------------------------------------------------------------------------
#
# Base path of the Spatial Registration Tool
export IP_PREFIX=$PWD
# ----------------------------------------------------------------------------------------
#
# Input Dataset - Scalar fields. VTK - "DATASET STRUCTURED_POINTS" or Raw-Format
export IP_DATA_IN=$IP_PREFIX/datasets/"Knochenprobe_2_77.vtk"
# export IP_DATA_IN=$IP_PREFIX/datasets/"Knochenprobe_2_7.vtk"
# ----------------------------------------------------------------------------------------
#
# Steering parameters
# ['0': Identity kernel, '1': Gaussian filter]
export IP_SELECT_K="1"  # Select Kernel
export IP_SIZE_K="3"    # Size of Kernel
export IP_GS="1.0"      # Sigma of Gauß Kernel
# ----------------------------------------------------------------------------------------
# Set a filename for the log file. Empty --> stdout / shell
export IP_LOG=$IP_PREFIX/$(date +"%m-%d-%y")"_log_ImageProcessing.txt"
# ----------------------------------------------------------------------------------------
#
# Feature deprecated with Version V0.2.0 due to performance issues. Now set in Source
# Set a Debugging level. "DEBUG" (1), "DEBUG_ALL" (2) or "PRODUCTION" (0). 
# export IP_DBG_LVL="2"
#
# Set Debug Flags for the mpi Compiler
# export IP_COMP_DBG="xterm -e gdb"
unset IP_COMP_DBG
# ----------------------------------------------------------------------------------------
#
# Set environment
# ----------------------------------------------------------------------------------------
if [ -n $IP_ARCH ]; then
    if [ $IP_ARCH = "julius" ]; then
        # --------------------------------------------------------------------------------
        # Make MPI available
        mpi_prefix="/opt/mpi/openmpi-4.1.0"
        #
        export PATH=${mpi_prefix}/bin:$PATH
        export LD_LIBRARY_PATH=${mpi_prefix}/lib:$LD_LIBRARY_PATH
        # --------------------------------------------------------------------------------
    elif [ $IP_ARCH = "vulcan "]; then
        # --------------------------------------------------------------------------------
        # Make MPI available
        module load mpi/openmpi/4.0.5-gnu-10.2.0
        #
        export PATH=${mpi_prefix}/bin:$PATH
        export LD_LIBRARY_PATH=${mpi_prefix}/lib:$LD_LIBRARY_PATH
        # --------------------------------------------------------------------------------
    elif [ $IP_ARCH == "hawk" ]; then
        # --------------------------------------------------------------------------------
        # Make MPI available
        module load mpt
    fi
else
    echo "Please define a valid architecture."
fi
# ----------------------------------------------------------------------------------------
