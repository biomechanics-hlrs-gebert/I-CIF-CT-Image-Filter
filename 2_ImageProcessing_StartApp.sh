#!/bin/bash
# -----------------------------------------------------------------------------
# Image Processing - Start Application
#
# Author:          Johannes Gebert «gebert@hlrs.de»
# Created:         25.04.2021
# Last edit:       09.05.2021
# ------------------------------------------------------------------------------
do_not_start=0
#
# Check whether all requirements are met
if [ -z $IP_PREFIX ]; then
  echo "Please source the environment file before:"
  echo "«source 1_ImageProcessing_Environment.sh»"
  do_not_start=1
fi
# -----------------------------------------------------------------------------
#
# Start program
if [ $do_not_start -eq 0 ]; then
    if [ -z $1 ]; then
		echo
        echo "Please specify the amount of processors required:"
        echo "./this_script «amount of processors»" && echo
		echo "At least 2, at max. the amount of processors of your computer."
		echo
		do_not_start=1
    else
		if [ $IP_ARCH == "julius" ]; then
			echo "Start Parallel Image Processing."
			echo
			mpirun $IP_COMP_DBG -np $1  $IP_PREFIX"/bin/HLRS_NUM_3D_Convolusional_Filtering_"$IP_ARCH"_"$IP_VRSN"_x86_64" 	 \
																										$IP_DATA_IN  \
																										$IP_LOG      \
																										$IP_MODE_K   \
																										$IP_SELECT_K \
																										$IP_SIZE_K   \
																										$IP_GS       \
																										$IP_VRSN	 \
																										$IP_CSV_TEX
																										# 
			status=$?			
		elif [ $IP_ARCH == "vulcan" ]; then
############################################################################
# VERY! raw. MODIFY BEFORE USE. 

#PBS -N run_Image_Processing
#PBS -l select=2:node_type=skl:mpiprocs=40
#PBS -l walltime=00:20:00

# Change to the direcotry that the job was submitted from
cd $PBS_O_WORKDIR

# using the INTEL MPI module
module load mpi/openmpi/4.0.5-gnu-10.2.0

# Launch the parallel mpi application (compiled with intel mpi) to the allocated compute nodes
mpirun -np 80 ./Image_Processing_julius_V0.1.0_x86_64 input.vtk out.log 1
############################################################################
		fi
    fi
fi
# ------------------------------------------------------------------------------
#
# Use this lines for basic postprocessing directives
if [ $do_not_start -eq 0 ]; then
	if [ $status -eq 0 ]; then
		echo
		echo "Parallel Image Processing finished successfully."
		echo "Postprocessing of Tex files started."
		echo

		which pdflatex > /dev/null 2> /dev/null
		if [ $? -eq 1 ]; then
			echo "Pdflatex not available."
			echo "Postprocessing stopped."
		else

			pdflatex -shell-escape -interaction nonstopmode -halt-on-error -file-line-error -output-directory=$PWD/tex/  $IP_CSV_TEX"_Filter_Histogram.tex" > $PWD/tex/"pdflatex_compile_"$ii"_.log"			
			
			if [ $? -eq 0 ]; then
				cd tex
				./TEX_clean_ignore.sh
				cd ..
			fi

			echo
			if [ $? -eq 1 ]; then
				echo "Postprocessing of Tex files crashed."
			else
				echo "Postprocessing of Tex files finished successfully."
			fi
		fi
	else
		echo
	    echo "Convolusional Filtering did not succeed. Postprocessing stopped."
	fi
fi

