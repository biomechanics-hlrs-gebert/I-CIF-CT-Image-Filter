#!/bin/bash
# -----------------------------------------------------------------------------
# Set the environment for the Julius system.
#    Julius I is a Whiskey Lake Intel(R) Core(TM) i5-8365U CPU @ 1.60GHz
#    laptop running Arch Linux x86_64, Kernel > 5.15.12-arch1-1, zsh > 5.8
#
# Author:    Johannes Gebert - HLRS - NUM - gebert@hlrs.de
# Created:   10.05.2021
# Last edit: 03.01.2022
# -----------------------------------------------------------------------------
# MPI environment
mpi_prefix=/opt/mpi/openmpi-I4-4.1.2
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
    metis_prefix=/opt/metis/metis-5.1.0
    export METIS_INCPATH=${metis_prefix}/include
    export METIS_LIBPATH=${metis_prefix}/lib
    #
    # ----------------------------------------
    # PETSc installation
    petsc_prefix=/opt/petsc/petsc-3.15
    export PETSC_INCPATH=${petsc_prefix}/include
    export PETSC_LIBPATH=${petsc_prefix}/lib
    export LD_LIBRARY_PATH=${petsc_prefix}/lib:$LD_LIBRARY_PATH
fi
# -----------------------------------------------------------------------------
# Define compile mode - Production or development
export compile_MODE=dev
# 
# -----------------------------------------------------------------------------
# Define std_out
export USE_STD_OUT=YES
#
# -----------------------------------------------------------------------------
# Gnu Debugger - make tmpi available / check prerequisites
tmpi_prefix="/opt/tmpi"
#
# -----------------------------------------------------------------------------
# Root is a git repo?
#
export PROVIDES_GIT="YES"
#
export PATH=${tmpi_prefix}:$PATH
#
tools=( gdb tmpi tmux mpirun )
#
dbg_err=0
#
if [ "$NO_OUTPUT" != "YES" ]; then
    for tool in "${tools[@]}"; do
        echo -n "-- "
        if ! which ${tool} ; then > /dev/null 2> /dev/null # (to suppress cmd line output)
            echo "-- Please provide ${yellow}${tool}${nc} to use gdb with mpi."
            dbg_err=1
        fi
    done
    #
    if [[ $dbg_err == 0 ]]; then
        echo "--"
        echo "-- Usage of the GNU Debugger:"
        echo "-- ${yellow}tmpi $1 gdb --args mpirun n_cpus binary-input-file${nc}"
        echo "-- After stopping gdb, [ctrl+b], [&], [y] and »exit« will get you "
        echo "-- back to the initial command line interface."
    fi
    #
    echo "-- "
fi