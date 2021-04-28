# --------------------------------------------------------------------------------------------------
# Makefile to build the Image Processing program
#
# Author:    Johannes Gebert «gebert@hlrs.de»
# Date:      25.04.2021
# Last edit: 25.04.2021
#
# For use of make visit: https://www.gnu.org/software/make/
# --------------------------------------------------------------------------------------------------
# Set Architecture to compile for
# hawk    - HPE Apollo
# vulcan  - NEC Cluster
# julius  - A Whiskey Lake Notebook, 4 cores, 16Gb memory, APU
trgt_arch = $(IP_ARCH)
trgt_vrsn = $(IP_VRSN)
# --------------------------------------------------------------------------------------------------
# Directories
mod_dir   = $(CURDIR)/mod/
obj_dir   = $(CURDIR)/obj/
lib_dir   = $(CURDIR)/lib/
bin_dir   = $(CURDIR)/bin/
f-src_dir = $(CURDIR)/f-src/
ext_f-src = $(CURDIR)/ext_f-src/
# --------------------------------------------------------------------------------------------------
# File extensions and suffixes
mod_ext   = .mod
obj_ext   = .o
sho_ext   = .so
f90_ext   = .f90
bin_suf   = x86_64
# --------------------------------------------------------------------------------------------------
clean_cmd = rm -f
# --------------------------------------------------------------------------------------------------
# Compilers
#ifeq($(strip $(trgt_arch)) ,"julius" )
  compiler = "mpif90"
#endif
#ifeq($(strip $(trgt_arch)) ,"hawk" )
#  compiler = "mpif90"
#endif
export compiler
# --------------------------------------------------------------------------------------------------
# Programming Environment - gnu, LLVM
PE = gnu
# --------------------------------------------------------------------------------------------------
# Compile flags GNU Compiler
ifeq ($(PE),gnu)
   c_flags_f90 = -J$(mod_dir) -I$(mod_dir) \
	             -g						   \
	             -O3					   \
	             -fbacktrace               \
                 -fbounds-check            \
	 			 -Wno-conversion           \
                 -Wall #                   \				# Diagnoses
#	             -fdefault-integer-8       \				# incompatible with ISO_FORTRAN_ENV
#	             -fdefault-real-8          \
#	             -finstrument-functions    \
#	             -fopenmp
endif
# --------------------------------------------------------------------------------------------------
# Executable
MAIN_bin = $(bin_dir)Image_Processing_$(trgt_arch)_$(trgt_vrsn)_$(bin_suf)
# --------------------------------------------------------------------------------------------------
# Generate objects
#
f-objects = $(obj_dir)mod_standards$(obj_ext)                          \
			$(obj_dir)mod_kernels$(obj_ext)                            \
			$(obj_dir)mod_stringmod$(obj_ext)                          \
			$(obj_dir)mod_file_routines$(obj_ext)  		               \
			$(obj_dir)mod_histogram$(obj_ext)  		               \
			$(obj_dir)main$(obj_ext)
# --------------------------------------------------------------------------------------------------
# Begin Building
all: $(MAIN_bin)

# -------------------------------------------------------------------------------------------------
# Standards Module
$(obj_dir)mod_standards$(obj_ext):$(f-src_dir)mod_standards$(f90_ext)
	@echo "---------------------------------------"
	@echo "-- Compiles: " $< "--------------------"
	$(compiler) $(c_flags_f90) -c $< -o $@
# -------------------------------------------------------------------------------------------------
# Module containing Convolutional matrices
$(obj_dir)mod_kernels$(obj_ext):$(f-src_dir)mod_kernels$(f90_ext)
	@echo "---------------------------------------"
	@echo "-- Compiles: " $< "--------------------"
	$(compiler) $(c_flags_f90) -c $< -o $@
# --------------------------------------------------------------------------------------------------
# Files routines module
$(obj_dir)mod_file_routines$(obj_ext):$(mod_dir)standards$(mod_ext) $(f-src_dir)mod_file_routines$(f90_ext)
	@echo "---------------------------------------"
	@echo "-- Compiles: " $(f-src_dir)mod_file_routines$(f90_ext) "--------------------"
	$(compiler) $(c_flags_f90) -c $(f-src_dir)mod_file_routines$(f90_ext) -o $@
# --------------------------------------------------------------------------------------------------
# Histogram module
$(obj_dir)mod_histogram$(obj_ext):$(mod_dir)standards$(mod_ext) $(mod_dir)file_routines$(mod_ext) \
																$(f-src_dir)mod_histogram$(f90_ext)
	@echo "---------------------------------------"
	@echo "-- Compiles: " $(f-src_dir)mod_histogram$(f90_ext) "--------------------"
	$(compiler) $(c_flags_f90) -c $(f-src_dir)mod_histogram$(f90_ext) -o $@
#  -------------------------------------------------------------------------------------------------
# External source to parse input
$(obj_dir)mod_stringmod$(obj_ext):$(mod_dir)standards$(mod_ext)	$(ext_f-src)stringmod$(f90_ext)
	@echo "---------------------------------------"
	@echo "-- Compiles: " $(ext_f-src)stringmod$(f90_ext) "--------------------"
	$(compiler) $(c_flags_f90) -c $(ext_f-src)stringmod$(f90_ext) -o $@
# --------------------------------------------------------------------------------------------------
# MAIN OBJECT
$(obj_dir)main$(obj_ext):$(mod_dir)standards$(mod_ext)			\
						 $(mod_dir)kernels$(mod_ext)			\
						 $(mod_dir)file_routines$(mod_ext)		\
 			             $(mod_dir)strings$(mod_ext)			\
			             $(mod_dir)histogram$(mod_ext)			\
						 $(f-src_dir)main$(f90_ext)
	@echo "---------------------------------------"
	@echo "-- Compiles: " $(f-src_dir)main$(f90_ext) "--------------------"
	$(compiler) $(c_flags_f90) -c $(f-src_dir)main$(f90_ext) -o $@
# ---------------------------------------------------------------------------------------------------
# Linking MAIN
$(MAIN_bin):$(f-objects)
	@echo "---------------------------------------"
	@echo "-- Linking MAIN binary"
	$(compiler) $(f-objects) -o $(MAIN_bin)

# ---------------------------------------------------------------------------------------------------
# Linking MAIN
	@echo "---------------------------------------"
	@echo "-- Successfully build all."
	@echo "---------------------------------------"


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
