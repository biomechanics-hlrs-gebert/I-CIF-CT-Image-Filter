#!/bin/bash
# ----------------------------------------------------------------------------------------
# Open-MPI download, build and install script.
#
# Author:          Johannes Gebert «gebert@hlrs.de»
# Created: on :    08.09.2020
# Last edit:       07.04.2021
# ----------------------------------------------------------------------------------------
#
# Requirements:
# GCC
# GNU Make
# tar
# wget
# Ensure, your OS is allowed to run this script:
# «sudo chmod +x Open-MPI_install.sh»
#
# ----------------------------------------------------------------------------------------
#
# usage: « ./<this script> »
# The Script will ask for sudo. It's required to create the target directory at /opt
# ----------------------------------------------------------------------------------------
#
# Modify only in case you know what you like to achieve :-)
# ----------------------------------------------------------------------------------------
#
# Specify the Version of MPI
VERSION="4.1.0"
# ----------------------------------------------------------------------------------------
#
# Check whether all programs and compilers are accesible
declare -a BIN=("gfortran" "gcc" "g++" "make" "tar" "wget")

ABORT=0

echo "----------------------------------------------------------------------------------------"
for i in "${BIN[@]}"
do
   which "$i" > /dev/null 2> /dev/null
   if [ $? -eq 1 ]; then
	echo "Please make «"$i"» available."
	ABORT=1
   fi
done

filename=openmpi-$VERSION.tar.gz

if [ $ABORT -eq 0 ]; then
	if [ ! -f $PWD/$filename ]; then
		wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-$VERSION.tar.gz
		if [ $? -eq 1 ]; then
			echo ""
			echo "Download failed. Please have a look for the tarball of Open-MPI and place it in the directory of this script."
			echo "Take a look at: https://www.open-mpi.org/software/ompi/v4.1/"
			ABORT=1
		fi
	fi
fi

if [ $ABORT -eq 0 ]; then

	PREPREF=/opt/mpi

	sudo mkdir -p $PREPREF >/dev/null 2>/dev/null &
	sudo chown $USER:$USER -R $PREPREF

	PREF=$PREPREF/openmpi-$VERSION

	tar -xvzf openmpi-$VERSION.tar.gz
	cd openmpi-$VERSION

	# Configure MPI.
	# --prefix=<installation directory, independent of this one.>
	# further options via ./configure --help
	export F77=gfortran
	export FC=gfortran
	export CC=gcc
	export CXX=g++
	export FFLAGS="-m64"
	export FCFLAGS="-m64"
	export CFLAGS="-m64"
	export CXXFLAGS="-m64"
	./configure --prefix=$PREF --enable-mpi-fortran=usempif08

	make && make install

	if [ $? -eq 0 ]; then
		rm -r ./openmpi-$VERSION* >/dev/null 2>/dev/null &
		echo ""
		echo "----------------------------------------------------------------------------------------"
		echo "Open-MPI successfully installed in:"
		echo ""
		echo $PREF
		echo ""
		echo "Append "$PREF"/bin to PATH and"
		echo "Append "$PREF"/lib to LD_LIBRARY_PATH."
		echo "Often accomplished via an environment script, specifically written to a program."
	else
		echo "The Installation of Open-MPI failed."
	fi
fi
echo "----------------------------------------------------------------------------------------"
