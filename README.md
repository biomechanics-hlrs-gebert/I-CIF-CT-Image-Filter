# Computed Tomography Image Filter
This program reads 3-Dimensional scalar fields out of STRUCTURED-POINTS \*.vtk files and filters them according to specific convolutional matrices/Kernels.  
It's tested with up to 160 Processors (4 nodes) on Vulcan. Turnaround time of roughly 200 Seconds while reading/computing/writing to storage on 18.5E09 Voxels of kind INTEGER2.

## Meta Template
Located in: 
```
./datasets/I-CTIF.meta.template
```
For use with previously used data sets:
```
cat ./datasets/I-CTIF.meta.template >> Your_Meta_File.meta
```
## [semantic versioning](https://semver.org):
Given a version number MAJOR.MINOR.PATCH, increment the:

* MAJOR version when you major Features (i.e. new way of image processing),
* MINOR version when you extend functionality (i.e. new kernels), and
* PATCH version when you make bug fixes.

## Requirements
* x86 64bit Hardware
* Linux x86 64Bit Installation with a Bash
* GNU Compiler Collection (GCC)
* An installation of Open-MPI
### Message Passing Interface 
Parallelization of the program is done with an API called MPI (Message Passing Interface).

Required: MPI - compiled with integer 4 and mpi_f08.

  1. [Open-mpi 4.1.0](https://www.open-mpi.org/software/ompi/v4.1/) on local systems. Other versions are not tested.
  2. [HPE-MPT on HLRS Hawk](https://kb.hlrs.de/platforms/index.php/MPI(Hawk))

The program may be ported to other architectures. Maybe not :-)

## Build
It's tested and therefore recommended to build and run the program as follows.
### Set up the Environment
```vim ./central_src/auxiliaries/system_environments/<system>.sh```
```source ./environment.sh <system>``` 

* Set an architecture/a system
  * Give the absolute base path of your mpi-installation
  * Alternatively give the proper module names of your compute cluster

### Run make:
Build the program:    ```make```
Create documentation: ```make docs```

### Uninstall:
```make clean && rm -r <your program directory>```

### Set up the parametrization
```./auxiliaries/pbs/CF.<system>.pbs```

## Usage
For example for testing on julius:
```
mpirun ./bin/ctif_v4.0.0_x86_64 -np 4  ./bin/ctif_v4.0.0_x86_64 ./datasets/SC00-0_tc_Dev_ctif_G3S11Sig10.meta
```
### Build the Diagrams after Filtering
```
cd ./datasets/
make tex=<*.tex file>
```

### Datasets
... are transfered via file exchange and are not pushed into the repository. 

### External Sources
Plain text headers are parsed via a [strings module](https://gbenthien.net/strings/index.html) by George Benthien from San Diego.
### Arbitrary
Use this program at your own risk.

