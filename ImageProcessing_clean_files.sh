#!/bin/bash
# -----------------------------------------------------------------------------
# Image Processing - Clean result files
#
# Author:          Johannes Gebert «gebert@hlrs.de»
# Created:         25.04.2021
# Last edit:       25.04.2021
# ------------------------------------------------------------------------------
#
# Add the files you like to delete by hand to get more control of compute intensive results
# ------------------------------------------------------------------------------
#
# Check whether trash is available
# Info: With trash, there is a fallback if you accidentally delete these files
which "trash" > /dev/null 2> /dev/null
if [ $? -eq 1 ]; then
    echo "Please make «trash» available."
    echo "You may download it here: https://github.com/andreafrancia/trash-cli"
    echo
    echo "Alternativly, use 'rm' to delete the files in this script."
else
    # If a file is not present, it will show an error message - which is fine.
    trash ./datasets/Knochenprobe_2_77_Kernel_0.vtk
    trash ./datasets/Knochenprobe_2_77_Kernel_1.vtk
    trash ./05-01-21_log_ImageProcessing.txt
fi
