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
# M-DDTC-Directly-Discretizing-Tensor-Computation / Struct process ...
# ... specific environments
#
# Basically a workaround to deal with centralized sources.
# Feasible due to the framework of the phd project.
# Simply checks for the struct process Fortran main file.
if [[ -f "f-src/struct_process.f90" ]]; then
    # ----------------------------------------
    # BLAS/LAPACK installation
    export LAPACK_LIBPATH=$PWD/lib/lapack
    #
    # ----------------------------------------
    # METIS installation
    module load numlib/metis/5.1.0-64bitint-gnu-10.3.0
    #
    metis_prefix=/opt/numlib/metis/5.1.0-64bitint-gnu-10.3.0
    export METIS_INCPATH=${metis_prefix}/include
    export METIS_LIBPATH=${metis_prefix}/lib
    #
    # ----------------------------------------
    # PETSc installation
    module load numlib/petsc/3.14.1-openmpi-4.0.5-gnu-10.2.0
    #
    petsc_prefix=/opt/numlib/petsc/3.14.1-openmpi-4.0.5-gnu-10.2.0
    export PETSC_INCPATH=${petsc_prefix}/include
    export PETSC_LIBPATH=${petsc_prefix}/lib
    export LD_LIBRARY_PATH=${petsc_prefix}/lib:$LD_LIBRARY_PATH
fi
# -----------------------------------------------------------------------------
#
# Define std_out
export USE_STD_OUT=YES
# -----------------------------------------------------------------------------
#
# Root is a git repo?
export PROVIDES_GIT=NO