
MAKEOVERRIDE=

#------------------------------------------------------------------------
# Directory names
#
ROOT=../crack
MODULAR=../Modular
DISL=../Disl

LAMMPS_SRC = ../lammps/src
LAMMPS_FORT2 = ../lammps/examples/COUPLE/fortran2
#LAMMPS_SRC = /home/rzamora/lammps-4May14/src
#LAMMPS_FORT2 = /home/rzamora/lammps-4May14/examples/COUPLE/fortran3

#------------------------------------------------------------------------
# Problem specific object files and their dependencies
#
OBJ1  = mesh.o 
FIL1  = mesh.f 
DEP1  = $(MODULAR)/mod_grain.f $(MODULAR)/mod_boundary.f \
		$(MODULAR)/mod_global.f $(OBJ3)
OBJ2  = field.o
FIL2  = field.f
DEP2  = $(MODULAR)/mod_global.f $(MODULAR)/mod_file.f \
		$(MODULAR)/mod_grain.f $(OBJ3)
OBJ3  = mod_crack.o
FIL3  = mod_crack.f


#----------------------------------------------------------------
# Fortran compile flags
# 
#OPTIM = -O1 #-g
OPTIM = -g
LIN_FF90="-lmpi -w -c $(OPTIM) -FI -pc80 -mp -prec-div -I$(MODULAR)"
LIN_FF90f="-lmpi -w -c $(OPTIM) -FR -pc80 -mp -prec-div -I$(MODULAR)"
#LIN_f90=ifort
LIN_f90=mpif90
LIN_LIBS=" -lm -L$(LAMMPS_SRC) -llammps_hive1 -L$(LAMMPS_FORT2) -llammps_fortran -L$(MODULAR) -lmdlr -lpthread -L$(DISL) -ldisl"
#LIN_LIBS="-lmpi -lm -L$(LAMMPS_SRC) -llammps_hive1 -L$(LAMMPS_FORT2) -llammps_fortran -L$(MODULAR) -lmdlr -lpthread -L$(DISL) -ldisl"

#F90=ifort
#----------------------------------------------------------------
# Do all of these
#

default:
	@echo "make linux"

linux:
	make  FF90=$(LIN_FF90) FF90f=$(LIN_FF90f) F90=$(LIN_f90) \
		LIBS=$(LIN_LIBS) crack 

linux_mod_crack:
	make FF90=$(LIN_FF90) FF90f=$(LIN_FF90f) F90=$(LIN_f90) \
                LIBS=$(LIN_LIBS) mod_crack.o



#----------------------------------------------------------------
# Targets 
#
crack: mesh.o field.o mod_crack.o $(MODULAR)/libmdlr.a $(DISL)/libdisl.a
	$(F90) -o crack $(MODULAR)/pmain.o mesh.o field.o mod_crack.o\
		$(LIBS)

mesh.o: $(DEP1) mesh.f
	$(F90) $(FF90) mesh.f
field.o: $(DEP2) field.f
	$(F90) $(FF90) field.f
mod_crack.o: mod_crack.f
	$(F90) $(FF90f) mod_crack.f


#------------------------------------------------------------------------
# Utilities
#
list:
	ls -l $(ROOT)/*.f $(MODULAR)/*.f

tar:
	tar cvf $(ROOT).tar $(ROOT)/*.f $(MODULAR)/*.f $(ROOT)/mak* $(MODULAR)/make

clean:
	-rm -f *.mod *.o bigcrack


