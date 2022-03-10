SYSTEM ?= zam044

## Check whether CUDA compiler is available
ENABLE_GPU := $(shell command -v nvcc 2> /dev/null)
#ENABLE_GPU := $() # Use this line to test non-CUDA version on a system with CUDA installed

# For reproducability compile the git key (hash) into the code
GITKEY=$(shell git log | head -1 | sed -e 's/commit //g')

# THISMACHINE := $(shell hostname)

EXECUTABLE = green
all: ${EXECUTABLE}

OBJECTS = \
	control.o \
	recorded_warnings.o \
	display_units.o \
	green.o


.PHONY: clean allclean
clean:
	-rm -f *.o *.so *~ *.d *.a

allclean: clean
	-rm -f ${EXECUTABLE}

CU     = nvcc
CC     = gcc
CXX    = g++

INC ?=
INC += -I../include
INC += -I../external
# INC += -I/usr/local/cuda/include


ERRFLAGS =
FEAFLAGS =
# FEAFLAGS += -D_Output_Units_Fixed
FEAFLAGS += -D_GIT_KEY=$(GITKEY)
FEAFLAGS += -D HAS_NO_CUDA
FEAFLAGS += -D HAS_RAPIDXML
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
OPTFLAGS += -g -pedantic -Wall
# gcc-options
# OPTFLAGS += -Wno-sign-compare -Wno-format -Wno-format-security
### -mtune=power8 -mcpu=power8 -mpower8-fusion -mpower8-vector -mvsx -maltivec -mdirect-move -fstrict-aliasing -fomit-frame-pointer -fno-schedule-insns


# LDFLAGS = -mkl
LDFLAGS = ## -lmkl_intel_thread

CUOPTFLAGS =
#### production mode
CUOPTFLAGS += -O3 --use_fast_math
#### debug mode
# CUOPTFLAGS += -g -O0 -G
# CUOPTFLAGS += -DDEBUGGPU

## verbose assembler
CUOPTFLAGS +=  -Xptxas -v

# for profiling
CUOPTFLAGS += -lineinfo
ifdef ENABLE_GPU
	FEAFLAGS += -DUSE_NVTX
	LDFLAGS  += -rdynamic -L/usr/local/cuda/lib64 -lcudart -lnvToolsExt
endif

EXTFLAGS ?=

ifeq ($(SYSTEM),zam044)
	EXTFLAGS += -D HAS_NO_MKL
	EXTFLAGS += -D HAS_FFTW
	LDFLAGS += -lm -lblas -llapack ### works on zam044 (MacBook Pro)
	LDFLAGS += -lfftw3
else
	## eq ($(THISMACHINE), "zam054")
	### works with g++ on zam054 (JSC workstation)
	LDFLAGS += -lm -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread 
endif

FLAGS = ${ERRFLAGS} ${FEAFLAGS} ${OPTFLAGS} ${EXTFLAGS}

GENCODE_SM60= -gencode arch=compute_60,code=sm_60 -gencode arch=compute_60,code=compute_60
GENCODE_SM70= -gencode arch=compute_70,code=sm_70 -gencode arch=compute_70,code=compute_70
GENCODE_FLAGS= ${GENCODE_SM60} ${GENCODE_SM70}

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
	${CXX} ${CFLAGS} ${INC} -MM $^ -MF $@

-include ${DEP}
