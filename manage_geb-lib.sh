#!/bin/bash
# ----------------------------------------------------------------------------------------
# Johannes Gebert - Doctoral project - initialize the subtree of the Geberts Library    
#
# Author:    Johannes Gebert - HLRS - NUM - gebert@hlrs.de
# Created:   27.12.2021
# Last edit: 27.02.2022
# ----------------------------------------------------------------------------------------
# Update the subtree
#
if  ! which git > /dev/null 2> /dev/null; then 
    echo "Git not found."
    echo "The program cannot get updates from this directory."
    echo "The program cannot compile if the »geb-lib« directory ist missing."
fi
#
if  ls -l "$PWD"/central_src > /dev/null 2> /dev/null; then 
    operation="pull"
else
    operation="add"
fi
#
git subtree $operation --prefix \
geb-lib git@github.com:biomechanics-hlrs-gebert/A-GLI-Geberts-Library.git \
main --squash

