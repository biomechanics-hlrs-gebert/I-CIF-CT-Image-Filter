# ------------------------------------------------------------------------------
# Makefile to build the centralized sources of the doctoral project of 
#
# Author:    Johannes Gebert - HLRS - NUM - gebert@hlrs.de
# Date:      27.12.2021
# Last edit: 27.12.2021
#
# For use of make visit: https://www.gnu.org/software/make/
# ------------------------------------------------------------------------------
# Handling of the subtree hash is not entirely clear yet (Dec 2021)
# 
# ifeq ($(PROVIDES_GIT),YES)
# # Get git hash https://jblevins.org/log/vc
# 	rev = $(shell git rev-parse HEAD)
# else
# 	rev = NO_GIT_REPOSITORY
# endif
# ------------------------------------------------------------------------------
long_name="centralized sources"
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
f90_compiler = "gfortran"
export f90_compiler
c_compiler = "gcc"
export c_compiler
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
# Generate objects
#
f-objects = $(obj_dir)mod_global_std$(obj_ext)\
			$(obj_dir)mod_strings$(obj_ext)\
			$(obj_dir)mod_messages_errors$(obj_ext) \
			$(obj_dir)mod_meta$(obj_ext)\
			$(obj_dir)mod_vtk_raw$(obj_ext)
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
# External source to parse input
$(obj_dir)mod_strings$(obj_ext):$(mod_dir)global_std$(mod_ext)	$(ext_f-src)strings$(f90_ext)
	@echo "----- Compiling " $(ext_f-src)strings$(f90_ext) " -----"
	$(f90_compiler) $(c_flags_f90) -c $(ext_f-src)strings$(f90_ext) -o $@
	@echo

# -----------------------------------------------------------------------------
#-- Error Handling Module 
$(obj_dir)mod_messages_errors$(obj_ext):$(mod_dir)global_std$(mod_ext) $(mod_dir)strings$(mod_ext) \
									$(f-src_dir)mod_messages_errors$(f90_ext)
	@echo "----- Compiling " $(f-src_dir)mod_messages_errors$(f90_ext) " -----"
	$(f90_compiler) $(c_flags_f90) -c $(f-src_dir)mod_messages_errors$(f90_ext) -o $@
	@echo 

# -----------------------------------------------------------------------------
#-- Meta Module 
$(obj_dir)mod_meta$(obj_ext):$(mod_dir)strings$(mod_ext) $(mod_dir)messages_errors$(mod_ext) \
							$(f-src_dir)mod_meta$(f90_ext)
	@echo "----- Compiling " $(f-src_dir)mod_meta$(f90_ext) " -----"
	$(f90_compiler) $(c_flags_f90) -c $(f-src_dir)mod_meta$(f90_ext) -o $@
	@echo 

# ------------------------------------------------------------------------------
# Files routines module
$(obj_dir)mod_vtk_raw$(obj_ext):$(mod_dir)global_std$(mod_ext)\
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