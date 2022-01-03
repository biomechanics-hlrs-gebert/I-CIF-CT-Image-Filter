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
    # module load scalapack
    export LAPACK_LIBPATH=$PWD/lib/lapack
    # 
    # export LAPACK_LIBPATH=/opt/hlrs/spack/rev-004_2020-06-17/scalapack/2.1.0-gcc-9.2.0-amna4d3j
    #
    # ----------------------------------------
    # METIS installation
    module load metis/5.1.0-int64
    #
    metis_prefix=/opt/hlrs/spack/rev-004_2020-06-17/metis/5.1.0-gcc-9.2.0-rdkkxlua
    export METIS_INCPATH=${metis_prefix}/include
    export METIS_LIBPATH=${metis_prefix}/lib
    #
    # ----------------------------------------
    # PETSc installation
    # module load petsc/3.12.2-int64
    #
    petsc_prefix=$PWD/petsc
    export PETSC_INCPATH=${petsc_prefix}/include
    export PETSC_LIBPATH=${petsc_prefix}/lib
    export LD_LIBRARY_PATH=${petsc_prefix}/lib:$LD_LIBRARY_PATH
fi
#
# ----------------------------------------
# Define std_out
export USE_STD_OUT=NO
#
# ----------------------------------------
# Root is a git repo?
export PROVIDES_GIT="NO"