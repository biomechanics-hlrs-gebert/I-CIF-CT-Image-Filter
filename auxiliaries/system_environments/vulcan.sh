#!/bin/bash
# -----------------------------------------------------------------------------
# Set the environment for the Vulcan system, a heterogenous cluster.
#
# Author:          Johannes Gebert - HLRS - NUM - gebert@hlrs.de
# Created:         25.12.2021
# Last edit:       25.12.2021
# -----------------------------------------------------------------------------
#
# MPI environment
module load mpi/openmpi/4.1.0-gnu-10.3.0 > /dev/null 2> /dev/null
#
mpi_prefix=/opt/hlrs/mpi/openmpi/4.0.5-gnu-10.2.0
export PATH=${mpi_prefix}/bin:$PATH
export LD_LIBRARY_PATH=${mpi_prefix}/lib:$LD_LIBRARY_PATH
# -----------------------------------------------------------------------------
#
# Define std_out
export USE_STD_OUT=NO
# -----------------------------------------------------------------------------
#
# Root is a git repo?
export PROVIDES_GIT=NO