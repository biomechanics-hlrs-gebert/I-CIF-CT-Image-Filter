#!/bin/bash
# ----------------------------------------------------------------------------------------
# Johannes Gebert - Doctoral project - initialize the subtree of the centralized sources
#
# Author:          Johannes Gebert - HLRS - NUM - gebert@hlrs.de
# Created:         27.12.2021
# Last edit:       28.12.2021
# ----------------------------------------------------------------------------------------
# Update the subtree
#
which git > /dev/null 2> /dev/null
if  [ $? -ne 0 ] ; then 
    echo "Git not found."
    echo "The program cannot get updates from this directory."
    echo "The program cannot compile if the »central_src« directory ist missing."
fi
#
ls -l $PWD/central_src > /dev/null 2> /dev/null
if  [ $? -ne 0 ] ; then 
    operation="add"
else
    operation="pull"
fi
#
git subtree $operation --prefix \
central_src git@github.com:biomechanics-hlrs-gebert/A-CESO-Central_Sources.git \
main --squash