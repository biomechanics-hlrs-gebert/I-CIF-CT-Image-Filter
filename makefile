# ------------------------------------------------------------------------------
# Makefile to build the CTIF-CT-Image-Filter
#
# Author:    Johannes Gebert - HLRS - NUM - gebert@hlrs.de
# Date:      25.04.2021
# Last edit: 25.12.2021
#
# For use of make visit: https://www.gnu.org/software/make/
# ------------------------------------------------------------------------------
trgt_vrsn="v4.0.0"
bin_name="ctif"
long_name="Computed Tomography Image Filter"
# ------------------------------------------------------------------------------
ifeq ($(PROVIDES_GIT),YES)
# Get git hash https://jblevins.org/log/vc
	rev = $(shell git rev-parse HEAD)
else
	rev = NO_GIT_REPOSITORY
endif
# -----------------------------------------------------------------------------
# Check for environment
check-env:
ifeq ($(IP_SYS),)
	@echo "-----------------------------------------------"
	@echo "-- Please source environment.sh <system> first."
	@echo "-----------------------------------------------"
else
	@echo "-----------------------------------------------"
	@echo "-- Environment to build for: "$(IP_SYS)
	@echo "-----------------------------------------------"
	$(MAKE) all
endif
# ------------------------------------------------------------------------------
# Build path
build_path = $(CURDIR)
export build_path
#
# ------------------------------------------------------------------------------
# Directories
mod_dir   = $(build_path)/mod/
obj_dir   = $(build_path)/obj/
lib_dir   = $(build_path)/lib/
bin_dir   = $(build_path)/bin/
f-src_dir = $(build_path)/f-src/
ext_f-src = $(build_path)/f-src/ext-src_
#
# Directory for documentation
doc_dir  = $(build_path)/doc/
html_dir = $(build_path)/html/
tex_dir  = $(build_path)/latex/
#
# ------------------------------------------------------------------------------
# File extensions and suffixes
mod_ext = .mod
obj_ext = .o
sho_ext = .so
f90_ext = .f90
bin_suf = _x86_64
# ------------------------------------------------------------------------------
clean_cmd = rm -f
# ------------------------------------------------------------------------------
# Compilers
#ifeq($(strip $(trgt_arch)) ,"julius" )
  compiler = "mpif90"
#endif
#ifeq($(strip $(trgt_arch)) ,"hawk" )
#  compiler = "mpif90"
#endif
export compiler
# ------------------------------------------------------------------------------
# Programming Environment - gnu, LLVM
PE = gnu
# ------------------------------------------------------------------------------
# Compile flags GNU Compiler
ifeq ($(PE),gnu)
   c_flags_f90 = -J$(mod_dir) -I$(mod_dir) \
				-g \
				-o \
				-O3 \
				-fbacktrace \
				-fbounds-check \
				-fbackslash \
				-Wno-conversion \
				-Wall
endif
# ------------------------------------------------------------------------------
# Executable
main_bin = $(bin_dir)$(bin_name)_$(trgt_vrsn)$(bin_suf)

# ------------------------------------------------------------------------------
# Generate objects
#
f-objects = $(obj_dir)mod_global_std$(obj_ext)\
			$(obj_dir)mod_kernels$(obj_ext)\
			$(obj_dir)mod_strings$(obj_ext)\
			$(obj_dir)mod_file_routines_mpi$(obj_ext)\
			$(obj_dir)mod_aux_routines_ip$(obj_ext)\
			$(obj_dir)ct_image_filter$(obj_ext)
# ------------------------------------------------------------------------------
# Begin Building
all: $(main_bin)

# ------------------------------------------------------------------------------
# Standards Module
$(obj_dir)mod_global_std$(obj_ext):$(f-src_dir)mod_global_std$(f90_ext)
	@echo "----- Compiling " $(f-src_dir)mod_global_std$(f90_ext) " -----"
	$(compiler) $(c_flags_f90) -c $(f-src_dir)mod_global_std$(f90_ext) -o $@
	@echo

# ------------------------------------------------------------------------------
# Module containing Convolutional matrices
$(obj_dir)mod_kernels$(obj_ext):$(f-src_dir)mod_kernels$(f90_ext)
	@echo "----- Compiling " $(f-src_dir)mod_kernels$(f90_ext)" -----"
	$(compiler) $(c_flags_f90) -c $(f-src_dir)mod_kernels$(f90_ext) -o $@
	@echo

# ------------------------------------------------------------------------------
# Files routines module
$(obj_dir)mod_file_routines_mpi$(obj_ext):$(mod_dir)global_std$(mod_ext) $(f-src_dir)mod_file_routines_mpi$(f90_ext)
	@echo "----- Compiling " $(f-src_dir)mod_file_routines_mpi$(f90_ext) " -----"
	$(compiler) $(c_flags_f90) -c $(f-src_dir)mod_file_routines_mpi$(f90_ext) -o $@
	@echo

# ------------------------------------------------------------------------------
# Histogram module
$(obj_dir)mod_aux_routines_ip$(obj_ext):$(mod_dir)global_std$(mod_ext) $(f-src_dir)mod_aux_routines_ip$(f90_ext)
	@echo "----- Compiling " $(f-src_dir)mod_aux_routines_ip$(f90_ext) " -----"
	$(compiler) $(c_flags_f90) -c $(f-src_dir)mod_aux_routines_ip$(f90_ext) -o $@
	@echo

# ------------------------------------------------------------------------------
# External source to parse input
$(obj_dir)mod_strings$(obj_ext):$(mod_dir)global_std$(mod_ext)	$(ext_f-src)strings$(f90_ext)
	@echo "----- Compiling " $(ext_f-src)strings$(f90_ext) " -----"
	$(compiler) $(c_flags_f90) -c $(ext_f-src)strings$(f90_ext) -o $@
	@echo

# ------------------------------------------------------------------------------
# Main object 
$(obj_dir)ct_image_filter$(obj_ext):$(mod_dir)global_std$(mod_ext)\
						 $(mod_dir)kernels$(mod_ext)\
						 $(mod_dir)file_routines_mpi$(mod_ext)\
 			             $(mod_dir)strings$(mod_ext)\
			             $(mod_dir)aux_routines_ip$(mod_ext)\
						 $(f-src_dir)ct_image_filter$(f90_ext)
	@echo "----- Compiling " $(f-src_dir)ct_image_filter$(f90_ext) " -----"
	$(compiler) $(c_flags_f90) -c $(f-src_dir)ct_image_filter$(f90_ext) -o $@

# -----------------------------------------------------------------------------
# Final Link step of MAIN -----------------------------------------------------
$(main_bin):$(f-objects)
	@echo "----------------------------------------------------------------------------------"
	@echo '--- Write revision and git info'
	@echo "CHARACTER(LEN = scl), PARAMETER :: revision = '$(trgt_vrsn)'" > $(f-src_dir)include_f90/revision_meta$(f90_ext)
	@echo "CHARACTER(LEN = scl), PARAMETER :: hash = '$(rev)'" >> $(f-src_dir)include_f90/revision_meta$(f90_ext)
	@echo "----------------------------------------------------------------------------------"
	@echo '--- Final link step of $(long_name) executable'
	@echo "----------------------------------------------------------------------------------"
	$(compiler) $(f-objects) -o $(main_bin)
	@echo
	@echo "----------------------------------------------------------------------------------"
	@echo "-- Successfully build all."
	@echo "----------------------------------------------------------------------------------"

help:
	@echo "----------------------------------------------------------------------------------"
	@echo "$(long_name) make targets"
	@echo "Regular:       »make (all)«   - Build the $(long_name)"
	@echo "Cleaning:      »make clean«   - Remove generated files, keep the config"
	@echo "Documentation: »make docs     - Build the html and the tex documentation"
	@echo "----------------------------------------------------------------------------------"

docs: 
	@echo "----------------------------------------------------------------------------------"
	@echo "-- Beginn buiding the documentation of the $(long_name)."
	@echo "----------------------------------------------------------------------------------"
	doxygen doc/doxy.conf
	$(MAKE) pdf -C $(tex_dir)  
	@echo "----------------------------------------------------------------------------------"
	@echo "-- Successfully build the documentation of the $(long_name)."
	@echo "----------------------------------------------------------------------------------"

cleandocs:
	@echo "----------------------------------------------------------------------------------"
	@echo "-- Cleaning html documentation"
	@echo "----------------------------------------------------------------------------------"
	$(clean_cmd) $(html_dir)/*
	@echo "----------------------------------------------------------------------------------"
	@echo "-- Cleaning tex documentation"
	@echo "----------------------------------------------------------------------------------"
	$(clean_cmd) $(tex_dir)/*
	@echo "----------------------------------------------------------------------------------"
	@echo "-- Documentation removed."
	@echo "----------------------------------------------------------------------------------"
	
clean:
	@echo "----------------------------------------------------------------------------------"
	@echo "-- Cleaning module directory"
	@echo "----------------------------------------------------------------------------------"
	$(clean_cmd) $(mod_dir)*$(mod_ext)
	@echo "----------------------------------------------------------------------------------"
	@echo "-- Cleaning object directory"
	@echo "----------------------------------------------------------------------------------"
	$(clean_cmd) $(f-objects)
	@echo "----------------------------------------------------------------------------------"
	@echo "-- Cleaning MAIN binary"
	@echo "----------------------------------------------------------------------------------"
	$(clean_cmd) $(MAIN_bin)
	@echo "----------------------------------------------------------------------------------"
	@echo "-- Cleaning completed."
	@echo "----------------------------------------------------------------------------------"
