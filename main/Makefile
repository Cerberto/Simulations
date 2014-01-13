################################################################################
#
#   Makefile to compile and link Fortran programs. (Hopefully) written for 
#	merging C routines.
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

###############################  main programs  ################################

MAIN =  test # lennard-jones
	#  test_vettori


###################  requided modules and external routines  ###################

EXTRAS = zeros

GLOBAL = kind  constants

ODE = numerov

SIMULATIONS = ljmod  ljmetropolis
	# metropolis
RANDOM =

CLUSTER = jackknife

COSEACASO = funzioni procedure_ex


MODULES =   $(GLOBAL) $(CLUSTER) $(RANDOM) # $(SIMULATIONS)
	# $(EXTRAS) $(ODE) $(COSEACASO)
	
	
######################  additional files for C routines  #######################

CFILES = ranlxs  ranlxd  gauss  ranluxf  start  utils  rintf


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

CINCLUDE = ../include  /usr/include

CC = gcc

CFLAGS = -Wall -O -std=gnu99 -g -DSSE2

CLIBS = m

CINCDIRS = $(addprefix -I,$(CINCLUDE))

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

FCFLAGS = -Wall -Wextra -C -g3 -fbacktrace -O2 -fbounds-check

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
	@ echo "  linking object files and producing executables"


# clean directory
rmpmk:
	@ -rm -rf *.o .tmp; \
		echo -e "  delete object and temporary files\n\n\t\tSTAY BOMBER!\n"

.PHONY: clean

clean: rmpmk
	@ rm -f *.mod *.MOD *~ $(PROGRAMS); \
		echo "  main directory cleaned"


################################################################################
