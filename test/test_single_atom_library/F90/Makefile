###
### Makefile to test the usage of libliveatom.a with Fortran
###

F90=gfortran
LD=g++
F90_FLAGS=
LD_FLAGS= 

### name of the execution file
EXECUTABLE= test_live_atom
OBJECTS=    test_live_atom.o

### directory containing source files
LIBDIR= ../../../src
LIBS= -L${LIBDIR} -lliveatom -lblas -llapack \
	-L/usr/local/Cellar/gcc/10.2.0_4/lib/gcc/10 -lgfortran

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
