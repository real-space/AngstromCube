
## Check whether CUDA compiler is available
ENABLEGPU := $(shell command -v nvcc 2> /dev/null)
#ENABLEGPU := $() # Use this line to test non-CUDA version on a system with CUDA installed

EXECUTABLE = aa3
all: ${EXECUTABLE}

OBJ = radial_grid.o radial_integrator.o radial_eigensolver.o element_configuration.o \
	radial_potential.o atom_core.o overlap.o single_atom.o angular_grid.o sho_radial.o main.o
ifdef ENABLEGPU
	OBJ += GPUdrivers.o
endif


.PHONY: clean allclean
clean:
	-rm -f *.o *~ *.d

allclean: clean
	-rm -f ${EXECUTABLE}

CU     = nvcc
### -G
CC     = gcc
CXX    = g++ ### icpc -mkl
##MPI
# CC  = mpicc
# CXX = mpic++

SRC ?= ./src
INC = -I/usr/local/cuda/include

EXTFLAGS ?=

ERRFLAGS = -Werror -Wall -Wfatal-errors
FEAFLAGS = ###-fopenmp -fno-omit-frame-pointer
# FEAFLAGS += -D NO_UNIT_TESTS ## for production

# normal mode
# OPTFLAGS = -Ofast -mtune=power8 -mcpu=power8 -mpower8-fusion -mpower8-vector -mvsx -maltivec -mdirect-move -fstrict-aliasing -fomit-frame-pointer -fno-schedule-insns
OPTFLAGS = -O2 ###-Ofast -fstrict-aliasing -fomit-frame-pointer -fno-schedule-insns

LDFLAGS =

CUOPTFLAGS =
#### production mode
CUOPTFLAGS += -O3 --use_fast_math
# CUOPTFLAGS += -D BENCH_ALL_VERSIONS
#### debug mode
# OPTFLAGS  = -O0 -g
# CUOPTFLAGS += -g -O0 -G
# CUOPTFLAGS += -DDEBUGGPU

## verbose assembler
CUOPTFLAGS +=  -Xptxas -v

# for profiling
CUOPTFLAGS += -lineinfo
ifdef ENABLEGPU
	FEAFLAGS += -DUSE_NVTX
	LDFLAGS  += -rdynamic -L/usr/local/cuda/lib64 -lcudart -lnvToolsExt
endif

FLAGS = ${ERRFLAGS} ${FEAFLAGS} ${OPTFLAGS} ${EXTFLAGS}

GENCODE_SM60= -gencode arch=compute_60,code=sm_60 -gencode arch=compute_60,code=compute_60
GENCODE_SM70= -gencode arch=compute_70,code=sm_70 -gencode arch=compute_70,code=compute_70
GENCODE_FLAGS= ${GENCODE_SM60} ${GENCODE_SM70}

CCFLAGS  = -std=gnu99 ${FLAGS}
# K80
#CUFLAGS  = -std=c++11 ${GENCODE_FLAGS} ${CUOPTFLAGS} $(patsubst %,-Xcompiler %,${FLAGS})
# P100
CUFLAGS  = -std=c++11  ${CUOPTFLAGS} $(patsubst %,-Xcompiler %,${FLAGS})
CUFLAGS  = -std=c++11  ${GENCODE_FLAGS} ${CUOPTFLAGS} $(patsubst %,-Xcompiler %,${FLAGS})
CXXFLAGS = -std=c++11 ${FLAGS}

# LDFLAGS += -lm -llapack 
LDFLAGS += -lm -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

DEP=$(OBJ:%.o=%.d)

${EXECUTABLE}: ${OBJ}
	${CXX} ${CXXFLAGS} $^ ${LDFLAGS}  -o $@

%.o:${SRC}/%.c %.d
	${CC} ${CCFLAGS} ${INC} -c $<
%.o:${SRC}/%.cu %.d
	${CU} ${CUFLAGS} ${INC} -c $<
%.o:${SRC}/%.cxx %.d
	${CXX} ${CXXFLAGS} ${INC} -c $<

%.o:%.c %.d
	${CC} ${CCFLAGS} ${INC} -c $<
%.o:%.cu %.d
	${CU} ${CUFLAGS} ${INC} -c $<
%.o:%.cxx %.d
	${CXX} ${CXXFLAGS} ${INC} -c $<

%.d:${SRC}/%.c
	${CC} ${CFLAGS} ${INC} -MM $^ -MF $@
%.d:${SRC}/%.cu
	${CC} -x c++ ${CFLAGS} ${INC} -MM $^ -MF $@
%.d:${SRC}/%.cxx
	${CXX} ${CXXFLAGS} ${INC} -MM $^ -MF $@

%.d:%.c
	${CC} ${CFLAGS} ${INC} -MM $^ -MF $@
%.d:%.cu
	${CC} -x c++ ${CFLAGS} ${INC} -MM $^ -MF $@
%.d:%.cxx
	${CXX} ${CXXFLAGS} ${INC} -MM $^ -MF $@

-include ${DEP}
