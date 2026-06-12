INCLUDE=.:/usr/local/include 

DIR_GEOMETRY = ../Geometry_Modules/Bulk

DIR_GENERAL = ../General_Routines
DIR_POLYMERS = ../Polymer_Modules
DIR_INTERACTIONS = ../Interaction_Modules

OBJS =	\
  mod_global.o \
  mod_mt.o \
  mod_fourier.o \
  mod_propagate.o \
  mod_iterate.o \
  sub_wrap.o \
  mod_FloryHuggins.o \
  mod_multiblock.o 

LIBS =-lfftw3	

CC = g++
CFLAGS = -O2
FC = gfortran
FFLAGS = -O
F90 = gfortran
F90FLAGS = -O3  -march=native -pipe -fomit-frame-pointer 
FCHECKFLAGS = -fbounds-check -pipe -fomit-frame-pointer 
LDFLAGS = -fopenmp -L/usr/lib:/usr/local/lib

all: 

main_multiblock: $(OBJS) main_multiblock.f90
	$(F90) -I$(INCLUDE) $(F90FLAGS) -c $@.f90
	$(F90) -I$(INCLUDE) $(LDFLAGS) -L/usr/lib:/usr/local/lib $(OBJS) $@.o $(LIBS) -o $@.exe


.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90 
	$(F90)  -I$(INCLUDE) $(F90FLAGS) -c $<

%.o: $(DIR_GENERAL)/%.f90
	$(F90)  -I$(INCLUDE) $(F90FLAGS) -c $<

%.o: $(DIR_GEOMETRY)/%.f90
	$(F90)  -I$(INCLUDE) $(F90FLAGS) -c $<

%.o: $(DIR_POLYMERS)/%.f90
	$(F90)  -I$(INCLUDE) $(F90FLAGS) -c $<

%.o: $(DIR_INTERACTIONS)/%.f90
	$(F90)  -I$(INCLUDE) $(F90FLAGS) -c $<

clean:
	rm -f *.o *.mod fort.*
