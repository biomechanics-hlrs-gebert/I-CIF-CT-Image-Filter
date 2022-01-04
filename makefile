# ------------------------------------------------------------------------------
# Makefile to build the centralized sources of the doctoral project of 
#
# Author:    Johannes Gebert - HLRS - NUM - gebert@hlrs.de
# Date:      27.12.2021
# Last edit: 29.12.2021
#
# For use of make visit: https://www.gnu.org/software/make/
# ------------------------------------------------------------------------------
long_name="Centralized Sources"
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
subtree_build_path = $(CURDIR)
export subtree_build_path
#
# ------------------------------------------------------------------------------
# Directories
mod_dir   = $(subtree_build_path)/mod/
obj_dir   = $(subtree_build_path)/obj/
c-src_dir = $(subtree_build_path)/c-src/
f-src_dir = $(subtree_build_path)/f-src/
ext_f-src = $(subtree_build_path)/f-src/ext-src_
#
# Directory for documentation
doc_dir  = $(subtree_build_path)/doc/
html_dir = $(subtree_build_path)/html/
tex_dir  = $(subtree_build_path)/latex/
#
# ------------------------------------------------------------------------------
# File extensions and suffixes
mod_ext = .mod
obj_ext = .o
sho_ext = .so
f90_ext = .f90
c_ext = .c
# ------------------------------------------------------------------------------
clean_cmd = rm -f
# ------------------------------------------------------------------------------
# Compilers
f90_compiler = "mpif90"
export f90_compiler
c_compiler = "gcc"
export c_compiler
# ------------------------------------------------------------------------------
# Programming Environment - gnu, LLVM
PE = gnu
# ------------------------------------------------------------------------------
# Compile mode - dev, prod - Defined in environment.sh
# compile_MODE = dev
# ------------------------------------------------------------------------------
# Compile flags GNU Compiler
# The subtree structure requires two directories containing modules. 
# In this case, the program root/mod directory addressed by the -J 
# http://www.hpc.icc.ru/documentation/intel/f_ug1/fced_mod.htm
ifeq ($(PE),gnu)
	f90_std_IJ     = -J$(mod_dir) -I$(subtree_build_path)
	f90_dev_flags  = 	-fdefault-integer-8 -fdefault-real-8 \
						-finstrument-functions -ggdb -o -O3 \
						-fbacktrace -fbounds-check \
						-Wno-conversion -Wall
	f90_prod_flags = 	-fdefault-integer-8 -fdefault-real-8 \
						-finstrument-functions -O3 -fbounds-check

	ifeq ($(compile_MODE),prod)
		c_flags_f90 = $(f90_std_IJ) $(f90_prod_flags)
	else
		c_flags_f90 =  $(f90_std_IJ) $(f90_dev_flags)
	endif
endif
# ------------------------------------------------------------------------------
# Generate objects
#
f-objects = $(obj_dir)mod_global_std$(obj_ext)\
			$(obj_dir)mod_strings$(obj_ext)\
			$(obj_dir)mod_math$(obj_ext)\
			$(obj_dir)mod_mechanical$(obj_ext)\
			$(obj_dir)mod_user_interaction$(obj_ext)\
			$(obj_dir)mod_meta$(obj_ext)\
			$(obj_dir)mod_vtk_raw$(obj_ext)\
			$(obj_dir)mod_formatted_plain$(obj_ext)

# ------------------------------------------------------------------------------
# Begin Building
all: $(f-objects)

# ------------------------------------------------------------------------------
# Standards Module
$(obj_dir)mod_global_std$(obj_ext):$(f-src_dir)mod_global_std$(f90_ext)
	@echo "----- Compiling " $(f-src_dir)mod_global_std$(f90_ext) " -----"
	$(f90_compiler) $(c_flags_f90) -c $(f-src_dir)mod_global_std$(f90_ext) -o $@
	@echo

# ------------------------------------------------------------------------------
# Mechanical Module
$(obj_dir)mod_mechanical$(obj_ext):$(mod_dir)global_std$(mod_ext)	$(f-src_dir)mod_mechanical$(f90_ext)
	@echo "----- Compiling " $(f-src_dir)mod_mechanical$(f90_ext) " -----"
	$(f90_compiler) $(c_flags_f90) -c $(f-src_dir)mod_mechanical$(f90_ext) -o $@
	@echo

# ------------------------------------------------------------------------------
# External source to parse input
$(obj_dir)mod_strings$(obj_ext):$(mod_dir)global_std$(mod_ext)	$(ext_f-src)strings$(f90_ext)
	@echo "----- Compiling " $(ext_f-src)strings$(f90_ext) " -----"
	$(f90_compiler) $(c_flags_f90) -c $(ext_f-src)strings$(f90_ext) -o $@
	@echo

# ------------------------------------------------------------------------------
# Math Module
$(obj_dir)mod_math$(obj_ext):$(mod_dir)global_std$(mod_ext)	\
							$(mod_dir)strings$(mod_ext)\
							$(f-src_dir)mod_math$(f90_ext)
	@echo "----- Compiling " $(f-src_dir)mod_math$(f90_ext) " -----"
	$(f90_compiler) $(c_flags_f90) -c $(f-src_dir)mod_math$(f90_ext) -o $@
	@echo
	
# -----------------------------------------------------------------------------
# Module for User Interaction
$(obj_dir)mod_user_interaction$(obj_ext):$(mod_dir)global_std$(mod_ext) $(mod_dir)strings$(mod_ext) \
									$(f-src_dir)mod_user_interaction$(f90_ext)
	@echo "----- Compiling " $(f-src_dir)mod_user_interaction$(f90_ext) " -----"
	$(f90_compiler) $(c_flags_f90) -c $(f-src_dir)mod_user_interaction$(f90_ext) -o $@
	@echo 

# -----------------------------------------------------------------------------
# Meta Module 
$(obj_dir)mod_meta$(obj_ext):$(mod_dir)strings$(mod_ext) $(mod_dir)user_interaction$(mod_ext) \
							$(f-src_dir)mod_meta$(f90_ext)
	@echo "----- Compiling " $(f-src_dir)mod_meta$(f90_ext) " -----"
	$(f90_compiler) $(c_flags_f90) -c $(f-src_dir)mod_meta$(f90_ext) -o $@
	@echo 

# -----------------------------------------------------------------------------
# Formatted Plain Module
$(obj_dir)mod_formatted_plain$(obj_ext):$(mod_dir)global_std$(mod_ext) $(mod_dir)math$(mod_ext)\
										$(f-src_dir)mod_formatted_plain$(f90_ext)
	@echo "----- Compiling " $(f-src_dir)mod_formatted_plain$(f90_ext) " -----"
	$(f90_compiler) $(c_flags_f90) -c $(f-src_dir)mod_formatted_plain$(f90_ext) -o $@
	@echo 

# ------------------------------------------------------------------------------
# Module vtk structured points and raw data
$(obj_dir)mod_vtk_raw$(obj_ext):$(mod_dir)global_std$(mod_ext) \
								$(mod_dir)user_interaction$(mod_ext) \
								$(f-src_dir)mod_vtk_raw$(f90_ext)
	@echo "----- Compiling " $(f-src_dir)mod_vtk_raw$(f90_ext) " -----"
	$(compiler) $(c_flags_f90) -c $(f-src_dir)mod_vtk_raw$(f90_ext) -o $@
	@echo

help:
	@echo "----------------------------------------------------------------------------------"
	@echo "Make targets"
	@echo "Regular:       »make (all)«   - Build all $(long_name)."
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
	@echo "-- Cleaning central_src module directory"
	@echo "----------------------------------------------------------------------------------"
	$(clean_cmd) $(mod_dir)*$(mod_ext)
	@echo "----------------------------------------------------------------------------------"
	@echo "-- Cleaning central_src object directory"
	@echo "----------------------------------------------------------------------------------"
	$(clean_cmd) $(obj_dir)*$(obj_ext)