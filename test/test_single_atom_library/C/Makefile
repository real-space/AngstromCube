###
### Makefile to test the usage of libliveatom.a with C
###

CC=gcc
LD=gcc
CC_FLAGS=
LD_FLAGS=

### name of the execution file
EXECUTABLE= test_live_atom
OBJECTS=    test_live_atom.o

### directory containing source files
LIBDIR= ../../../src
LIBS= -L${LIBDIR} -lliveatom -lblas -llapack

### directory containing header files
INCDIR= ../../../include

###############################
## RULES ######################
###############################
# the pattern how to create .o-object out of a .CC-source file

%.o: %.c
	$(CC) -I$(INCDIR) $(CC_FLAGS) -c $<

###############################
### TARGETS ###################
###############################
### first target (is taken when calling make without a target argument)

$(EXECUTABLE): $(OBJECTS)
	$(LD) $(LD_FLAGS) $^ $(LIBS)  -o $(EXECUTABLE)
