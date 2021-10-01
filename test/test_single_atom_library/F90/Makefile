###
### Makefile to test the usage of libliveatom.a with Fortran
###

F90=gfortran
# LD=g++
LD=gfortran
F90_FLAGS=
LD_FLAGS= 

### name of the execution file
EXECUTABLE= test_live_atom
OBJECTS=    test_live_atom.o

### directory containing source files
LIBDIR= ../../../src
### if the linker is c++ and we link libliveatom.a
LIBS= -L${LIBDIR} -lliveatom -lblas -llapack -L/usr/local/Cellar/gcc/10.2.0_4/lib/gcc/10 -lgfortran
### if the linker is fortran and we link libliveatom.so
LIBS= -L${LIBDIR} -lliveatom

### directory containing header files
INCDIR= ../../../include

###############################
## RULES ######################
###############################
# the pattern how to create .o-object out of a .F90-source file

%.o: %.F90
	$(F90) -I$(INCDIR) $(F90_FLAGS) -c $<

###############################
### TARGETS ###################
###############################
### first target (is taken when calling make without a target argument)

$(EXECUTABLE): $(OBJECTS)
	$(LD) $(LD_FLAGS) $^ $(LIBS)  -o $(EXECUTABLE)
