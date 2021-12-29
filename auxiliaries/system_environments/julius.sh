#!/bin/bash
# -----------------------------------------------------------------------------
# Set the environment for the Julius system.
#    Julius I is a Whiskey Lake Intel(R) Core(TM) i5-8365U CPU @ 1.60GHz
#    laptop running Arch Linux x86_64, Kernel > 5.14.12-arch1-1, zsh > 5.8
#
# Author:    Johannes Gebert - HLRS - NUM - gebert@hlrs.de
# Created:   10.05.2021
# Last edit: 29.12.2021
# -----------------------------------------------------------------------------
#
# MPI environment
mpi_prefix=/opt/mpi/openmpi-NO_F08-4.1.0
export PATH=${mpi_prefix}/bin:$PATH
export LD_LIBRARY_PATH=${mpi_prefix}/lib:$LD_LIBRARY_PATH
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
