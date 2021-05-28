#!/bin/bash
# ----------------------------------------------------------------------------------------
# Set Environment for Image Processing
#
# Author:          Johannes Gebert - HLRS - NUM «gebert@hlrs.de»
# Created:         10.05.2021
# Last edit:       12.05.2021
# ----------------------------------------------------------------------------------------
#
# Set current Version
export IP_VRSN="V3.2.1"
# ----------------------------------------------------------------------------------------
#
#
#
# ----------------------------------------------------------------------------------------
#
usage ()
{
    echo 
    echo "Usage: source source2set_Environment.sh «environment»"
    echo 
    echo "Environments available:"
    echo "hawk    - HPE Apollo"
    echo "vulcan  - NEC Cluster"
    echo "julius  - A Whiskey Lake Notebook, 4 cores, 16Gb memory, APU"
    echo
}
#
manual ()
{
    echo
    echo "This Script sets the environment to compile or run the 3D Convolusional Filter."
    echo "Sometimes it's referred as «Image Processing». Currently, the terms are interchangeable."
    echo
    echo "You may set a Version number within this script, which has to stay consistent to Git."
    echo "Current Version: "$IP_VRSN
    usage
}
#
err ()
{
    echo 
    echo "Environment for Architecture "$IP_ARCH" was not sourced successfully." 
    echo "Please check the system architecture or check for other errors."
    usage  
    success=0
}
#
#------------------------------------------------------------------------------
#
export IP_PREFIX=$PWD
#
if [ -z $1 ]; then
    manual
    # exit 1
else 
    export IP_ARCH=$1
    #
    success=1
    #
    # Set environment
    case $IP_ARCH in

    julius)
        mpi_prefix="/opt/mpi/openmpi-NO_F08-4.1.0"
        export PATH=${mpi_prefix}/bin:$PATH
        export LD_LIBRARY_PATH=${mpi_prefix}/lib:$LD_LIBRARY_PATH
        ;;

    vulcan)
        module load mpi/openmpi/4.1.0-gnu-10.3.0 > /dev/null 2> /dev/null

        if [ $? -eq 0 ]; then
            export PATH=${mpi_prefix}/bin:$PATH
            export LD_LIBRARY_PATH=${mpi_prefix}/lib:$LD_LIBRARY_PATH
        else
            err
        fi
        ;;

    hawk)
        module load mpt > /dev/null 2> /dev/null

        if [ $? -ne 0 ]; then
            err
        fi
        ;;

    *)
        err
        ;;
    esac
fi
#
if [ $success -eq 1 ]; then
    echo "Environment for Architecture "$IP_ARCH" sourced."
fi
# ----------------------------------------------------------------------------------------