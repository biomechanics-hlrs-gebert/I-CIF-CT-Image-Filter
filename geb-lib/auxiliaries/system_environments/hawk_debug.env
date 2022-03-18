#!/bin/bash
# -----------------------------------------------------------------------------
# Set the environment for the Hawk system, a heterogenous cluster.
#
# Author:    Johannes Gebert - HLRS - NUM - gebert@hlrs.de
# Created:   25.12.2021
# Last edit: 28.01.2022
# -----------------------------------------------------------------------------
#
# Allocate compute node(s) before sourcing this script
# qsub -I -l select=1:node_type=rome:ncpus=128:mpiprocs=128 -l walltime=07:59:59
#
# Load VNC
/sw/vulcan-CentOS7/hlrs/tools/VirtualGL/2.6.5/bin/vncserver -rfbauth $HOME/.vnc/passwd
#
# Export Display.
export DISPLAY=:1
echo "-- Check whether DISPLAY=:1 is correct."
echo "-- After a second log in, this might be :2."
#
# Load debugger
module load forge
#
# Start debugger
# forge >/dev/null 2>/dev/null &
#
# On client (r35c3t4n2:1 --> compute_node:display)
# ./vncviewer -via hpcgeber@hawk-login04.hww.hlrs.de r35c3t4n2:1
# -----------------------------------------------------------------------------
#
# Load MPI
module load openmpi/4.0.5
#
mpi_prefix=/opt/hlrs/non-spack/mpi/openmpi/4.0.5-gcc-9.2.0/
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
    #
    metis_prefix=$PWD/lib/metis/metis-5.1.0
    export METIS_INCPATH=${metis_prefix}/include
    export METIS_LIBPATH=${metis_prefix}/lib
    #
    # ----------------------------------------
    # PETSc installation
    #
    petsc_prefix=$PWD/lib/petsc/petsc-3.15
    export PETSC_INCPATH=${petsc_prefix}/include
    export PETSC_LIBPATH=${petsc_prefix}/lib
    export LD_LIBRARY_PATH=${petsc_prefix}/lib:$LD_LIBRARY_PATH
fi
#
# ----------------------------------------
# Define std_out
export USE_STD_OUT=YES
#
# ----------------------------------------
# Root is a git repo?
export PROVIDES_GIT=NO
#
# ----------------------------------------
# Load stats
watch -n 2 "echo 'Connect to:'; echo vncviewer -via hpcgeber@hawk-login04.hww.hlrs.de $(uname -a | cut -d ' ' -f 2)$DISPLAY; echo; echo; free -h --mega; echo; echo; qstat -a; echo; echo; ls -l datasets"
