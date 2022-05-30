
## Check whether CUDA compiler is available
ENABLE_GPU := $(shell command -v nvcc 2> /dev/null)
#ENABLE_GPU := $() # Use this line to test non-CUDA version on a system with CUDA installed

# For reproducability compile the git key (hash) into the code
GITKEY=$(shell git log | head -1 | sed -e 's/commit //g')

HOSTNAME := $(shell hostname)

EXECUTABLE = green
all: ${EXECUTABLE}

OBJECTS = \
	control.o \
	recorded_warnings.o \
	display_units.o \
	green.o


.PHONY: clean allclean
clean:
	-rm -f *.o *~ *.d

allclean: clean
	-rm -f ${EXECUTABLE}

CU     = nvcc
CC     = gcc
CXX    = g++

INC ?=
INC += -I../include
INC += -I../external

ERRFLAGS =
FEAFLAGS =
# FEAFLAGS += -D_Output_Units_Fixed
FEAFLAGS += -D_GIT_KEY=$(GITKEY)
FEAFLAGS += -D HAS_RAPIDXML

## use the tfQMRgpu library to invert the shifted Hamiltonian
FEAFLAGS += -D HAS_TFQMRGPU
INC += -I../../tfQMRgpu

# ERRFLAGS = -Werror -Wall -Wfatal-errors
# FEAFLAGS += -fopenmp -fno-omit-frame-pointer

## for production: eliminate the unit tests
# FEAFLAGS += -D NO_UNIT_TESTS

## include also code which is not meant for release,
## use github.com/BR903/cppp -UDEVEL to do partial preprocessing
FEAFLAGS += -D DEVEL

# normal mode
# OPTFLAGS = -Ofast
# OPTFLAGS = -O2
OPTFLAGS = -O0
OPTFLAGS += -g -Wall
# OPTFLAGS += -pedantic
# gcc-options
# OPTFLAGS += -Wno-sign-compare -Wno-format -Wno-format-security

# LDFLAGS = -mkl
LDFLAGS = ## -lmkl_intel_thread

CUOPTFLAGS =
#### production mode
CUOPTFLAGS += -O3 --use_fast_math
#### debug mode
# CUOPTFLAGS += -g -O0 -G
# CUOPTFLAGS += -DDEBUGGPU

## verbose GPU assembler
# CUOPTFLAGS += -Xptxas -v
CUOPTFLAGS += --expt-relaxed-constexpr

# for profiling
CUOPTFLAGS += -lineinfo
ifdef ENABLE_GPU
	INC += -I/usr/local/cuda/include
	FEAFLAGS += -DUSE_NVTX
	LDFLAGS  += -L/usr/local/cuda/lib64 -lcudart -lnvToolsExt -lcurand  #-rdynamic
else
	FEAFLAGS += -D HAS_NO_CUDA
endif

EXTFLAGS ?=

FLAGS = ${ERRFLAGS} ${FEAFLAGS} ${OPTFLAGS} ${EXTFLAGS}

GENCODE60= -gencode arch=compute_60,code=sm_60 -gencode arch=compute_60,code=compute_60 ## for P100
GENCODE70= -gencode arch=compute_70,code=sm_70 -gencode arch=compute_70,code=compute_70 ## for V100
GENCODE80= -gencode arch=compute_80,code=sm_80 -gencode arch=compute_80,code=compute_80 ## for A100
GENCODE_FLAGS=
# GENCODE_FLAGS += ${GENCODE60}
GENCODE_FLAGS += ${GENCODE70}
# GENCODE_FLAGS += ${GENCODE80}

CCFLAGS  = -std=gnu99 ${FLAGS}
# K80
#CUFLAGS  = -std=c++11 ${GENCODE_FLAGS} ${CUOPTFLAGS} $(patsubst %,-Xcompiler %,${FLAGS})
# P100
CUFLAGS  = -std=c++11 ${CUOPTFLAGS} $(patsubst %,-Xcompiler %,${FLAGS})
CUFLAGS  = -std=c++11 ${GENCODE_FLAGS} ${CUOPTFLAGS} $(patsubst %,-Xcompiler %,${FLAGS})
CXXFLAGS = -std=c++11 ${FLAGS}

ifdef ENABLE_GPU
else
	CU = $(CXX)
	CUFLAGS = $(CXXFLAGS)
endif


DEP=$(OBJECTS:%.o=%.d)

${EXECUTABLE}: ${OBJECTS}
	${CU} ${CUFLAGS} $^ ${LDFLAGS}  -o $@

%.o:%.c %.d
	${CC} ${CCFLAGS} ${INC} -c $<
%.o:%.cxx %.d
	${CXX} ${CXXFLAGS} ${INC} -c $<
%.o:%.cu %.d
	${CU} ${CUFLAGS} ${INC} -c $<

%.d:%.c
	${CC} ${CFLAGS} ${INC} -MM $^ -MF $@
%.d:%.cxx
	${CXX} ${CXXFLAGS} ${INC} -MM $^ -MF $@
%.d:%.cu
	${CU} ${CUFLAGS} ${INC} -MM $^ -MF $@

-include ${DEP}
