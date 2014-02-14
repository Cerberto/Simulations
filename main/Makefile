################################################################################
#
#   Makefile to compile and link Fortran and C programs.
#
#  "make" compiles and links the specified main programs and modules (in C)
#  using the specified libraries (if any), and produces the executables;
#  by default it deletes also object files after producing executables.
#
#  "make fp" does the same as "make" for Fortran programs.
# 
#  "make clean" removes all files apart from source codes.
#
################################################################################

all: rmxeq cmpcf mkcxeq rmpmk 
.PHONY: all

SHELL=/bin/bash


#################################  C programs  #################################

MAIN = linpot_lattice  varmin


RANDOM = ranlxs  ranlxd  gauss  ranluxf

SIMULATIONS = metropolis

CLUSTER = jackknife

START = start  utils

EXTRAS = rintf  zeros

MODULES =  $(RANDOM) $(SIMULATIONS) $(CLUSTER) $(START) $(EXTRAS)


#############################  Fortran programs  ###############################

FMAIN = lennard-jones

FGLOBAL = kind  constants

FSIMULATIONS = ljmod  ljmetropolis

FCLUSTER = jackknifef

FMODULES = $(FGLOBAL) $(FCLUSTER) $(FSIMULATIONS)


##############################  set search paths  ##############################

# search path for modules
MDIR = ../modules
VPATH = $(MDIR)/global:$(MDIR)/cluster:$(MDIR)/simulations:$(MDIR)/start:\
	$(MDIR)/random:$(MDIR)/extras


###############################  for C programs  ###############################

INCPATH = ../include  /usr/include

LIBPATH = /usr/local/lib /usr/lib

CC = gcc

FLAGS = -Wall -O -std=gnu99 -pedantic -g -DSSE2

LIBS = m

PROGRAMS = $(MODULES) $(MAIN)

OBJECTS = $(addsuffix .o,$(MODULES))

INCDIRS = $(addprefix -I,$(INCPATH))

LDFLAGS = $(addprefix -L,$(LIBPATH)) $(addprefix -l,$(LIBS))

-include $(addsuffix .d,$(PROGRAMS))


# rule to make dependencies
$(addsuffix .d,$(PROGRAMS)): %.d: %.c Makefile
	@ $(CC) -MM -ansi $(INCDIRS) $< -o $@


# rule to compile source programs
$(addsuffix .o,$(PROGRAMS)): %.o: %.c Makefile
	@ $(CC) $< -c $(FLAGS) $(INCDIRS) -o $@


$(MAIN): %: %.o $(OBJECTS) Makefile
	@ $(CC) $< $(OBJECTS) $(CFLAGS) $(LDFLAGS) -o $@


############################  for Fortran programs  ############################

FINCPATH = /usr/include ../include

FLIBPATH = /usr/local/lib

FC = gfortran

FCFLAGS = -C -g3 -fbacktrace -O3 -fbounds-check -Wall -Wextra

FLIBS = lapack  blas

FPROGRAMS = $(FMODULES) $(FMAIN)

FINCDIRS = $(addprefix -I,$(FINCPATH))

FOBJECTS = $(addsuffix .o,$(FMODULES)) $(addsuffix .o,$(MODULES))

FLDFLAGS = $(addprefix -L,$(FLIBPATH)) $(addprefix -l,$(FLIBS))


# compile source programs
$(addsuffix .o,$(FPROGRAMS)): %.o: %.f90 Makefile
	@ $(FC) $(FLDFLAGS) $(FCFLAGS) $(FINCDIRS) -c $<

# link object files and produce executables
$(FMAIN): %: %.o $(FOBJECTS)
	@ $(FC) $(FCFLAGS) $(FLDFLAGS) -o $@ $^


############################  List of instructions  ############################

# remove old executables and old error log file
rmxeq:
	@ -rm -f $(MAIN) $(FMAIN); \
        echo "  delete old executables"


cmpcf: $(addsuffix .o,$(PROGRAMS))

cmpff: $(addsuffix .o,$(FPROGRAMS))

mkcxeq: $(MAIN)
	@ echo "  link object files and produce executables (C source codes)"

mkfxeq: $(FMAIN)
	@ echo "  link object files and produce executables (Fortran source codes)"

# remove files except for executables and .mod
rmpmk:
	@ -rm -rf *.d *.o .tmp; \
		echo -e "  delete object and temporary files\n\n\t\tSTAY BOMBER!\n"

# clean the main directory
.PHONY: clean
clean: rmpmk
	@ rm -f *.mod *.MOD *~ $(PROGRAMS) $(FPROGRAMS); \
		echo "  main directory cleaned"

# compile only Fortran programs
.PHONY: fp
fp: rmxeq cmpff mkfxeq rmpmk

################################################################################
