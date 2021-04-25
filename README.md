# Image Processing - Convolusional Filtering of Scalar Fields.

![Architecture](https://img.shields.io/badge/Architecture-x86-green)
![OS](https://img.shields.io/badge/Linux-64Bit-green)
![version](https://img.shields.io/badge/version-0.1.0-red)
![Contributors](https://img.shields.io/badge/HLRS-NUM-blue)
![Contributors](https://img.shields.io/badge/Contributors-2-blue)


This program reads 3-Dimensional scalar fields out of STRUCTURED-POINTS *.vtk files and filters them according to specific convolutional matrices/Kernels.

Currently, it's under development. Terms like "Scalar Fields", "Voxel Grids" and "Array (of Field Quantities)" are context dependent but in some ways interchangable synonyms.

- [x] All documentation is done in :us: to quickly accomodate non-native speakers.
- [x] Updates done to get to publication-readiness.
- [ ] Publish repository for demonstrating current progress to partners and HLRS employees.

## Datasets and other large files...
... are transfered via file exchange and are not pushed into the repository. 

## Primary target of this program
Simply put, Image Processing aims at filtering Computed Tomography Scans (basically scalar fields) according to our need of exposing a proper histogram for thresholding and binarizing in subsequent programs.

## Development of new features
All new features or parts of the tool are developed via branches. At any time, a version which compiles and executes and which does not necessarily crash is kept on the Repo's master branch.

Developing new features and calulations with C or Fortran may render a cumbersome undertaking. Utilizing an interpreted language or interpreting framework like MatLab, Octave, Python or sth. else is welcomed. However, please push *all* of the development files into this repository/testing/A... directory structure.

### We apply *a slightly modified* [semantic versioning](https://semver.org):

Given a version number MAJOR.MINOR.PATCH, increment the:

* MAJOR version when you make incompatible changes,
* MINOR version when you add functionality in a backwards compatible manner, and
* PATCH version when you make backwards compatible bug fixes.
## Requirements

* x86 64bit Hardware
* Linux x86 64Bit Installation with a BASh
* GNU Compiler Collection (GCC) including gcc/gfortran
* An installation of Open MPI - run the script in the project's root directory.
### Optional
[DrawIO](https://sourceforge.net/projects/drawio-desktop.mirror/) for documentary purposes.

DrawIO does not require an installation. To access it via Shell without using an absolute path (example):

```
sudo mkdir /opt/drawio
cp ~/Downloads/drawio-x86_64-14.5.1.AppImage !$
sudo ln -s /opt/drawio/drawio-x86_64-14.5.1.AppImage /usr/bin/drawio
```
### Message Passing Interface 
Parallelization of the program is done with an API called MPI (Message Passing Interface). Currently, Open-MPI leads this development. However, adoptions to HPE-MPT for use on hawk may be required at some stage of the project.

Required: MPI - compiled with integer 4 and mpi_f08 - simply run ```./Open-MPI_install.sh```

  1. [Open-mpi >=4.1.0](https://www.open-mpi.org/software/ompi/v4.1/) on local systems. Other versions are not tested.
  2. [HPE-MPT on HLRS Hawk](https://kb.hlrs.de/platforms/index.php/MPI(Hawk))

The program may be ported to other architectures. Maybe not :-)

## Build/Install
It's tested and therefore recommended to build and rund the program as follows.
### Set up the Environment variables within 
```./1_ImageProcessing_Environment.sh```

1. Set an Architecture
2. Define the Hardware Architecture; optional:

   1. Give the absolute base path of your mpi-installation
   2. Alternatively give the proper module names of your compute cluster

### Run make:
```
source 1_ImageProcessing_Environment.sh
make
```
### Uninstall:
```
make clean
rm -r «your program directory»
```

## Usage
### Set up the parametrization with proper parameters within 
```source 1_ImageProcessing_Environment.sh```

1. (Set an Architecture)
2. Set a base path - usually the one this README.md is located in
3. Give the absolute path of an input dataset
5. Define steering parameters
6. Set a filename to write the log to
7. Define a Debugging level 

   * 0 = Production
   * 1 = Standard development debugging level
   * 2 = Detailed output of some subroutines and steps.

```./2_ImageProcessing_StartApp.sh```

You may need to manually remove result files. It's strongly recommended not to delete these data automatically, due to results of lots of compute hours, which may be lost.

It is strongly recommended to use an amount of processors by a power of 2. In example 2, 4, 8, 16, ...
The Minimum amount of processors is 2.
To get the best result on HLRS Hawk, you may want to exploit spatial locality via topology aware scheduling. Therefore, an amount of processors by a factor of 8 may be best. In example 8, 64, ...

However this is an assumption, which is not tested by end of April 2021.

### Create Testcases
Create new testcases in ./testing.
Some of the Testprograms are written in Fortran and compiled with a BASh-Script (not with a makefile).

```
sudo chmod +x «compile_testcase.f90»
./«compile_testcase.f90»
```
## Developers
| Surname   | Name     | Department   |
|-----------| ---------| -------------|
| Schnabel  | Benjamin | HLRS - NUM   |
| Gebert    | Johannes | HLRS - NUM   |
## Limits
None specifically detected.
## Arbitrary
Use this program at your own risk.