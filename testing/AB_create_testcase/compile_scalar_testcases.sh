#!/bin/bash
gfortran  -g -O3 -fbacktrace -fbounds-check -o "create_scalar_testcases_x86_64" create_vtk_testcases.f90
