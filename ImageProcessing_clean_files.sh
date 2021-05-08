#!/bin/bash
# -----------------------------------------------------------------------------
# Image Processing - Clean result files
#
# Author:          Johannes Gebert «gebert@hlrs.de»
# Created:         25.04.2021
# Last edit:       08.05.2021
# ------------------------------------------------------------------------------
#
# Add the files you like to delete by hand to get more control of compute intensive results
# ------------------------------------------------------------------------------
#
# Choose the program to delete files with.
program="trash"
# program="rm"
# ------------------------------------------------------------------------------
#
# List the files you like to delete with paths relative to this directory.
declare -a files=(  "./datasets/Knochenprobe_2_77_Kernel_0.vtk"   \
                    "./datasets/Knochenprobe_2_77_Kernel_1.vtk"   \
                    "./datasets/Knochenprobe_2_77_Kernel_3.vtk"   \
                    "./datasets/hk1_B60f_cCT_Kernel_3.vtk"        \
                    "$(date +"%m-%d-%y")_log_ImageProcessing.txt" \
                    "./tex/$(date +"%m-%d-%y")_ImageProcessing_hist_PRE__FILTER.csv" \
                    "./tex/$(date +"%m-%d-%y")_ImageProcessing_hist_POST_FILTER.csv" \
                    "./tex/$(date +"%m-%d-%y")_ImageProcessing_Filter_Histogram.tex" \
                    "./tex/$(date +"%m-%d-%y")_ImageProcessing_Filter_Histogram.pdf" \
                    )
# ------------------------------------------------------------------------------
if [[ -f "./tex/TEX_clean_ignore.sh" ]]; then
    cd tex
    bash ./TEX_clean_ignore.sh
    cd ..
fi
# ------------------------------------------------------------------------------
which $program > /dev/null 2> /dev/null
if [ $? -eq 1 ]; then
    echo "Please make «$program» available."
    if [ "$program" == "trash" ]; then
        echo "You may download it here: https://github.com/andreafrancia/trash-cli"
    fi
else
    for i in "${files[@]}"
    do
        $program "$i" > /dev/null 2> /dev/null
        if [[ -f "$i" ]]; then
            if [ $? -eq 1 ]; then
                    echo $i" was not deleted as expected."
                    ABORT=1
            fi
#        else
#            echo $i" does not exist." 
        fi
    done
fi
# ------------------------------------------------------------------------------
