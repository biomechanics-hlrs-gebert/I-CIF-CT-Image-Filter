#!/bin/bash
# -----------------------------------------------------------------------------
# Set the environment for the Hawk system, a heterogenous cluster.
#
# Author:          Johannes Gebert - HLRS - NUM - gebert@hlrs.de
# Created:         25.12.2021
# Last edit:       25.12.2021
# -----------------------------------------------------------------------------
#
# Unload MPT ------------------------------
module unload mpt/2.23
#
# Unload MPT ------------------------------
module load gcc/9.2.0
#
# MPI environment ------------------------
module load openmpi/4.0.4
#
mpi_prefix=/opt/hlrs/non-spack/mpi/openmpi/4.0.4-gcc-9.2.0/
export PATH=${mpi_prefix}/bin:$PATH
export LD_LIBRARY_PATH=${mpi_prefix}/lib:$LD_LIBRARY_PATH
#
# ----------------------------------------
# Define std_out
export USE_STD_OUT=NO
#
# ----------------------------------------
# Root is a git repo?
export PROVIDES_GIT="NO"