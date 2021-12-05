# 3D Convolusional Filtering

![Architecture](https://img.shields.io/badge/Architecture-x86-green)
![OS](https://img.shields.io/badge/Linux-64Bit-green)
![version](https://img.shields.io/badge/version-3.4.1-green)
![Contributors](https://img.shields.io/badge/HLRS-NUM-blue)

This program reads 3-Dimensional scalar fields out of STRUCTURED-POINTS \*.vtk files and filters them according to specific convolutional matrices/Kernels.  
It's tested with up to 160 Processors (4 nodes) on Vulcan. Turnaround time of roughly 200 Seconds while reading/computing/writing to storage on 18.5E09 Voxels of kind INTEGER2.

## Table of contents
- [Intent](#primary)
- [Development](#development)
- [Requirements](#requirements)
- [Build](#build)
- [Usage](#usage)
- [Additional Information](additional-information)

## Intent of this program
Simply put, 3D Convolusional Filtering aims at Computed Tomography Scans (basically scalar fields) for our need of exposing a proper histogram for thresholding and binarizing in subsequent programs.

## Development
Developing new features and calculations with C or Fortran may render a cumbersome undertaking. Utilizing an interpreted language or interpreting framework like MatLab, Octave, Python or sth. else is welcomed. However, please push *all* of the development files into this repository/testing/A... directory structure.

Furthermore, on julius, an environment and *.run script was adopted to run the program via GDB. It requires tmux and tmpi as it synchronizes the commands of all gdb instances.
Many thanks to Arno Mayrhofer and his contributors: [tmpi](https://github.com/Azrael3000/tmpi) 
### Testcases
Create new testcases in ./testing.
Some of the Testprograms are written in Fortran and compiled with a BASh-Script (not with a makefile).

```
sudo chmod +x compile_testcase.f90
./compile_testcase.f90
```
### We apply *a slightly modified* [semantic versioning](https://semver.org):

Given a version number MAJOR.MINOR.PATCH, increment the:

* MAJOR version when you major Features (i.e. new way of image processing),
* MINOR version when you extend functionality (i.e. new kernels), and
* PATCH version when you make bug fixes.
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

  1. [Open-mpi 4.1.0](https://www.open-mpi.org/software/ompi/v4.1/) on local systems. Other versions are not tested.
  2. [HPE-MPT on HLRS Hawk](https://kb.hlrs.de/platforms/index.php/MPI(Hawk))

The program may be ported to other architectures. Maybe not :-)

## Build
It's tested and therefore recommended to build and run the program as follows.
### Set up the Environment variables within 
```./source2set_Environment.sh```

1. Set an Architecture
2. Define the Hardware Architecture; optional:

   1. Give the absolute base path of your mpi-installation
   2. Alternatively give the proper module names of your compute cluster

### Run make:
```
source ./source2set_Environment.sh »Architecture«
make
```
### Uninstall:
```make clean && rm -r »your program directory«```

### NEC aurora TSUBASA
Build on an interactive of NEC aurora TSUBASA. Mpinfort, which is necessary to compile for this architecture, is available on the aurora nodes. These nodes must be used interactively to build the binary. If the binary already was compiled, a regular qsub to the vector queue suffices.

```
qsub -q vector -l select=1:node_type=aurora:mpiprocs=24,walltime=3600 -I
source source2set_Environment.sh aurora
make
```

## Usage
It's recommended to use the BASh scripts to control program flow. However, manual control of the program is possible. To run the binary, you have to source the Environment file, too.

### Set up the parametrization with proper parameters within 
```./bin/CF.vulcan.pbs```

1. Give the absolute path of an input dataset
2. Define steering parameters
3. Define the compute nodes

### Start the program

```./bin/HLRS_NUM_3D_Convolusional_Filtering.run```

You may need to manually remove result files. It's strongly recommended not to delete these data automatically, due to results of lots of compute hours, which may be lost.

It is strongly recommended to use an amount of processors by a power of 2. In example 2, 4, 8, 16, ...
The Minimum amount of processors is 2.
To get the best result on HLRS Hawk, you may want to exploit spatial locality via topology aware scheduling. Therefore, an amount of processors by a factor of 8 may be best. In example 8, 64, ...

However this is an assumption, which is not tested by end of April 2021.

### Datasets
... are transfered via file exchange and are not pushed into the repository. 


## Additional Information
### Developers
<table>
<thead>
  <tr>
    <th>Surname</th>
    <th>Name</th>
    <th>Department</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td>Schnabel</td>
    <td>Benjamin</td>
    <td>HLRS - NUM</td>
  </tr>
  <tr>
    <td>Gebert</td>
    <td>Johannes</td>
    <td>HLRS - NUM</td>
  </tr>
</tbody>
</table>


### External Sources
Plain text headers are parsed via a [strings module](https://gbenthien.net/strings/index.html) by George Benthien from San Diego.
### Arbitrary
Use this program at your own risk.

