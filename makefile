# ------------------------------------------------------------------------------
# Makefile to build the CTIF-CT-Image-Filter
#
# Author:    Johannes Gebert - HLRS - NUM - gebert@hlrs.de
# Date:      25.04.2021
# Last edit: 30.12.2021
#
# For use of make visit: https://www.gnu.org/software/make/
# ------------------------------------------------------------------------------
bin_name="ctif"
long_name="Computed Tomography Image Filter"
# ------------------------------------------------------------------------------
ifeq ($(PROVIDES_GIT),YES)
# Get git hash https://jblevins.org/log/vc
# rev = $(shell git describe --tags --always)
	rev = $(shell git rev-parse HEAD)
	trgt_vrsn = $(shell git describe --tags --abbrev=0)
else
	rev = NO_GIT_REPOSITORY
	trgt_vrsn = ""
endif
# -----------------------------------------------------------------------------
# Check for environment
check-env:
ifeq ($(SYS_ENV),)
	@echo "-----------------------------------------------"
	@echo "-- Please source environment.sh <system> first."
	@echo "-----------------------------------------------"
else
	@echo "-----------------------------------------------"
	@echo "-- Environment to build for: "$(SYS_ENV)
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
# st: "Subtree" - A git procedure to inherit another repository as some sort of
# submodule. https://gist.github.com/SKempin/b7857a6ff6bddb05717cc17a44091202
st_path= $(build_path)/central_src/
#
st_obj_dir = $(st_path)/obj/
st_mod_dir = $(st_path)/mod/
st_f-src_dir = $(st_path)/f-src/
#
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
# Compile mode - dev, prod
compile_MODE = dev
# ------------------------------------------------------------------------------
# Compile flags GNU Compiler
# The subtree structure requires two directories containing modules. 
# In this case, the program root/mod directory addressed by the -J 
# http://www.hpc.icc.ru/documentation/intel/f_ug1/fced_mod.htm
ifeq ($(PE),gnu)
	f90_std_IJ     = -J$(mod_dir) -I$(st_mod_dir)
	f90_dev_flags  = -ggdb -o -O3 -fbacktrace -fbounds-check -fbackslash -Wno-conversion -Wall
	f90_prod_flags = -O3 -fbounds-check

	ifeq ($(compile_MODE),prod)
		c_flags_f90 = $(f90_std_IJ) $(f90_prod_flags)
	else
		c_flags_f90 =  $(f90_std_IJ) $(f90_dev_flags)
	endif
endif
# ------------------------------------------------------------------------------
# Executable
main_bin = $(bin_dir)$(bin_name)_$(trgt_vrsn)$(bin_suf)

# ------------------------------------------------------------------------------
# Generate objects
#
f-objects = $(st_obj_dir)mod_global_std$(obj_ext)\
			$(st_obj_dir)mod_strings$(obj_ext)\
			$(st_obj_dir)mod_math$(obj_ext)\
			$(st_obj_dir)mod_user_interaction$(obj_ext) \
			$(st_obj_dir)mod_meta$(obj_ext) \
			$(st_obj_dir)mod_vtk_raw$(obj_ext)\
			$(st_obj_dir)mod_formatted_plain$(obj_ext) \
			$(obj_dir)mod_kernels$(obj_ext)\
			$(obj_dir)mod_histogram_routines$(obj_ext)\
			$(obj_dir)ct_image_filter$(obj_ext)

# ------------------------------------------------------------------------------
# Build the st directory first
st: 
	$(MAKE) all -C $(st_path)
	@echo 

# ------------------------------------------------------------------------------
# Begin Building
all: st $(main_bin)  

# ------------------------------------------------------------------------------
# Histogram module
$(obj_dir)mod_histogram_routines$(obj_ext):$(st_mod_dir)global_std$(mod_ext) \
											$(st_mod_dir)meta$(mod_ext) \
											$(st_mod_dir)formatted_plain$(mod_ext) \
											$(f-src_dir)mod_histogram_routines$(f90_ext)
	@echo "----- Compiling " $(f-src_dir)mod_histogram_routines$(f90_ext) " -----"
	$(compiler) $(c_flags_f90) -c $(f-src_dir)mod_histogram_routines$(f90_ext) -o $@
	@echo


# ------------------------------------------------------------------------------
# Module containing convolutional kernels
$(obj_dir)mod_kernels$(obj_ext):$(st_mod_dir)math$(mod_ext) $(f-src_dir)mod_kernels$(f90_ext)
	@echo "----- Compiling " $(f-src_dir)mod_kernels$(f90_ext)" -----"
	$(compiler) $(c_flags_f90) -c $(f-src_dir)mod_kernels$(f90_ext) -o $@
	@echo


# ------------------------------------------------------------------------------
# Main object 
$(obj_dir)ct_image_filter$(obj_ext):$(st_mod_dir)global_std$(mod_ext)\
									$(mod_dir)kernels$(mod_ext)\
									$(st_mod_dir)vtk_meta_data$(mod_ext)\
 			            			$(st_mod_dir)strings$(mod_ext)\
			            			$(mod_dir)histogram_routines$(mod_ext)\
									$(f-src_dir)ct_image_filter$(f90_ext)
	@echo "----- Compiling " $(f-src_dir)ct_image_filter$(f90_ext) " -----"
	$(compiler) $(c_flags_f90) -c $(f-src_dir)ct_image_filter$(f90_ext) -o $@

# --------------------------------------------------------------------------------------------------
# Export revision
export_revision:
	@echo "----------------------------------------------------------------------------------"
	@echo '-- Write revision and git info'
	@echo "CHARACTER(LEN=scl), PARAMETER :: longname = '$(long_name)'" > $(st_f-src_dir)include_f90/revision_meta$(f90_ext)
	@echo "CHARACTER(LEN=scl), PARAMETER :: revision = '$(trgt_vrsn)'" >> $(st_f-src_dir)include_f90/revision_meta$(f90_ext)
	@echo "CHARACTER(LEN=scl), PARAMETER :: hash = '$(rev)'" >> $(st_f-src_dir)include_f90/revision_meta$(f90_ext)
	@echo "----------------------------------------------------------------------------------"

# -----------------------------------------------------------------------------
# Final Link step of MAIN
$(main_bin): export_revision $(f-objects)
	@echo "----------------------------------------------------------------------------------"
	@echo '-- Final link step of $(long_name) executable'
	@echo "----------------------------------------------------------------------------------"
	$(compiler) $(f-objects) -o $(main_bin)
	@echo
	@echo "----------------------------------------------------------------------------------"
	@echo "-- Successfully build all."
	@echo "----------------------------------------------------------------------------------"

help:
	@echo "----------------------------------------------------------------------------------"
	@echo "-- $(long_name) make targets"
	@echo "-- Regular:  »make (all)«    - Build the $(long_name)"
	@echo "-- Cleaning: »make clean«    - Remove build files, keep the central_src"
	@echo "-- Cleaning: »make cleanall« - Remove all build files."
	@echo "-- Docs:     »make docs      - Build the html and the tex documentation."
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
	$(clean_cmd) $(obj_dir)*$(obj_ext)
	@echo "----------------------------------------------------------------------------------"
	@echo "-- Cleaning MAIN binary"
	@echo "----------------------------------------------------------------------------------"
	$(clean_cmd) $(main_bin)
	
cleanall: clean
	@echo "----------------------------------------------------------------------------------"
	@echo "-- Cleaning central_src st"
	@echo "----------------------------------------------------------------------------------"
	$(MAKE) clean -C $(st_path)