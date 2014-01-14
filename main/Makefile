################################################################################
#
#   Makefile to compile and link Fortran and C programs.
#
#  "make" compiles and links the specified main programs and modules
#  using the specified libraries (if any), and produces the executables;
#  by default it deletes also object files after producing executables.
# 
#  "make clean" removes all files apart from source codes.
#
################################################################################

all: rmxeq cmpcf cmpsc mkxeq rmpmk 
.PHONY: all

SHELL=/bin/bash

#############################  Fortran programs  ###############################

MAIN = lennard-jones


# requided modules and external routines

GLOBAL = kind  constants

ODE = numerov

SIMULATIONS = ljmod  ljmetropolis  ljroutines

RANDOM = 

CLUSTER = jackknife

MODULES =   $(GLOBAL) $(CLUSTER) $(RANDOM) $(SIMULATIONS)
	
	
#################################  C programs  #################################

CMAIN =

CRANDOM = ranlxs  ranlxd  gauss  ranluxf

CSTART = start  utils

CEXTRA = rintf

CFILES =  $(CRANDOM) $(CSTART) $(CEXTRA)


##############################  set search paths  ##############################

# search path for modules
MDIR = ../modules
VPATH = $(MDIR)/global:$(MDIR)/cluster:$(MDIR)/simulations:$(MDIR)/start:\
	$(MDIR)/random:$(MDIR)/extras

# additional include directories
INCPATH = /usr/include

# additional libraries to be included 
LIBPATH = /usr/local/lib


###############################  for C programs  ###############################

CINCPATH = ../include  /usr/include

CC = gcc

CFLAGS = -Wall -O -std=gnu99 -g -DSSE2

CLIBS = m

CINCDIRS = $(addprefix -I,$(CINCPATH))

CLDFLAGS = $(addprefix -L,$(LIBPATH)) $(addprefix -l,$(CLIBS))

-include $(addsuffix .d,$(CFILES))


# rule to make dependencies
$(addsuffix .d,$(CFILES)): %.d: %.c Makefile
	@ $(CC) -MM -ansi $(CINCDIRS) $< -o $@


# rule to compile source programs
$(addsuffix .o,$(CFILES)): %.o: %.c Makefile
	@ $(CC) $< -c $(CFLAGS) $(CINCDIRS) -o $@


############################  for Fortran programs  ############################

FC = gfortran

FCFLAGS = -C -g3 -fbacktrace -O2 -fbounds-check # -Wall -Wextra

LIBS = lapack  blas

PROGRAMS = $(MODULES) $(MAIN)

INCDIRS = $(addprefix -I,$(INCPATH))

OBJECTS = $(addsuffix .o,$(MODULES)) $(addsuffix .o,$(CFILES))

LDFLAGS = $(addprefix -L,$(LIBPATH)) $(addprefix -l,$(LIBS))


# compile source programs
$(addsuffix .o,$(PROGRAMS)): %.o: %.f90 Makefile
	@ $(FC) $(LDFLAGS) $(FCFLAGS) $(INCDIRS) -c $<

# link object files and produce executables
$(MAIN): %: %.o $(OBJECTS)
	@ $(FC) $(FCFLAGS) $(LDFLAGS) -o $@ $^


############################  List of instructions  ############################

# remove old executables and old error log file
rmxeq:
	@ -rm -f $(MAIN); \
        echo "  delete old executables"


# produce executables
cmpcf: $(addsuffix .o,$(CFILES))

cmpsc: $(addsuffix .o,$(PROGRAMS))

mkxeq: $(MAIN)
	@ echo "  link object files and produce executables"


# clean directory
rmpmk:
	@ -rm -rf *.d *.o .tmp; \
		echo -e "  delete object and temporary files\n\n\t\tSTAY BOMBER!\n"

.PHONY: clean

clean: rmpmk
	@ rm -f *.mod *.MOD *~ $(PROGRAMS); \
		echo "  main directory cleaned"


################################################################################
