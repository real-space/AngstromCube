# This file is part of AngstromCube under MIT License

### This makefile should only be used if CMake fails!

HOSTNAME := $(shell hostname)

# NO_UNIT_TESTS=1

# For reproducability compile the git key (hash) into the code
GITKEY=$(shell git log | head -1 | sed -e 's/commit //g')

EXECUTABLE = a43
all: ${EXECUTABLE}

OBJDIR=obj_a43

OBJATOM = $(addprefix $(OBJDIR)/, \
	control.o \
	recorded_warnings.o \
	radial_grid.o \
	radial_integrator.o \
	radial_potential.o \
	radial_eigensolver.o \
	exchange_correlation.o \
	scattering_test.o \
	chemical_symbol.o \
	display_units.o \
	sigma_config.o \
	atom_core.o \
	sho_unitary.o \
	angular_grid.o \
	single_atom.o )

OBJ = $(addprefix $(OBJDIR)/, \
	action_plan.o \
	dyadic_plan.o \
	kinetic_plan.o \
	geometry_input.o \
	geometry_analysis.o \
	self_consistency.o \
	iterative_poisson.o \
	parallel_poisson.o \
	bessel_transform.o \
	dense_solver.o \
	multi_grid.o \
	mpi_parallel.o \
	potential_generator.o \
	parallel_potential.o \
	poisson_solver.o \
	plane_wave.o \
	sho_basis.o \
	sho_overlap.o \
	sho_hamiltonian.o \
	energy_contour.o \
	sho_projection.o \
	sho_potential.o \
	structure_solver.o \
	green_parallel.o \
	green_action.o \
	green_input.o \
	green_kinetic.o \
	green_dyadic.o \
	green_experiments.o \
	green_sparse.o \
	green_solver.o \
	green_function.o \
	green_memory.o \
	green_potential.o \
	green_tests.o \
	load_balancer.o \
	a43_main.o)

.PHONY: clean allclean
clean:
	-rm -f *.o *.so *~ *.d *.a

allclean: clean
	-rm -f ${EXECUTABLE}

CC     = gcc
CXX   ?= g++
##MPI
# CC  = mpicc
# CXX = mpic++

INC ?=
INC += -I../include
INC += -I../external/rapidxml
INC += -I../external/rapidjson

ERRFLAGS =
FEAFLAGS =
FEAFLAGS += -D _GIT_KEY=$(GITKEY)
# FEAFLAGS += -D NO_UNIT_TESTS
# FEAFLAGS += -D _Output_Units_Fixed

# NDEBUG will disable all the assert() checks from <cassert> 
# FEAFLAGS += -D NDEBUG

FEAFLAGS += -D HAS_NO_MPI
FEAFLAGS += -D HAS_NO_OMP
FEAFLAGS += -D HAS_NO_CUDA
FEAFLAGS += -D HAS_BITMAP_EXPORT
FEAFLAGS += -D HAS_SINGLE_ATOM
# FEAFLAGS += -D USE_SIMPLE_RANGE_TRUNCATION
# FEAFLAGS += -D USE_RECIPROCAL_RADIAL_GRID
# FEAFLAGS += -D LARGE_ANGULAR_GRIDS
# FEAFLAGS += -D GENERAL_CELL
FEAFLAGS += -D HAS_TFQMRGPU
FEAFLAGS += -D HAS_RAPIDXML
FEAFLAGS += -D HAS_RAPIDJSON
FEAFLAGS += -D HAS_LAPACK
# ERRFLAGS = -Werror -Wall -Wfatal-errors
# FEAFLAGS += -fopenmp -fno-omit-frame-pointer

ifdef NO_UNIT_TESTS
	## for production: eliminate the unit tests
	FEAFLAGS += -D NO_UNIT_TESTS
endif

## include also code which is not meant for release,
## use github.com/BR903/cppp -UDEVEL to do partial preprocessing
FEAFLAGS += -D DEVEL

# normal mode
OPTFLAGS=
# OPTFLAGS += -Ofast
# OPTFLAGS += -O2
OPTFLAGS += -O0
OPTFLAGS += -g -pedantic -Wall
# gcc-options
OPTFLAGS += -Wno-format
# OPTFLAGS += -Wno-format-security
OPTFLAGS += -Wno-sign-compare

# LDFLAGS = -mkl
LDFLAGS = ## -lmkl_intel_thread


EXTFLAGS ?=

ifeq ($(HOSTNAME),MacBook-Pro-4)
	EXTFLAGS += -D HAS_NO_MKL
	LDFLAGS += -lm
	LDFLAGS += -lblas
	# OPTFLAGS += -fopenmp ### works only with gcc
	LDFLAGS += -llapack
	EXTFLAGS += -D HAS_FFTW
  # INC += -I/opt/homebrew/Cellar/fftw/3.3.10_1/include  ### needed by gcc
	LDFLAGS += -lfftw3 -L/opt/homebrew/Cellar/fftw/3.3.10_1/lib
else
	# on JUSUF:   module --force purge && module load Stages/2022 GCC/11.2.0 OpenMPI/4.1.2 imkl/2021.4.0
	LDFLAGS += -lm
	LDFLAGS += -lmkl_intel_lp64
	LDFLAGS += -lmkl_intel_thread
	LDFLAGS += -lmkl_core
	LDFLAGS += -liomp5
	LDFLAGS += -lpthread
endif

FLAGS = ${ERRFLAGS} ${FEAFLAGS} ${OPTFLAGS} ${EXTFLAGS}

CXXFLAGS = -std=c++11 ${FLAGS}


DEP=$(OBJ:%.o=%.d) $(OBJATOM:%.o=%.d)

${EXECUTABLE}: ${OBJATOM} ${OBJ}
		${CXX} ${CXXFLAGS} $^ ${LDFLAGS}  -o $@

libliveatom.a: ${OBJATOM} $(OBJDIR)/library_kind_static.o
		ar rcsv $@ $^

libliveatom.so: ${OBJATOM} $(OBJDIR)/library_kind_dynamic.o
		${CXX} ${CXXFLAGS} -fPIC -shared $^ ${LDFLAGS} -o $@

$(OBJDIR)/%.o:%.cxx $(OBJDIR)/%.d
		${CXX} ${CXXFLAGS} ${INC} -c $< -o $@

$(OBJDIR)/%.d:%.cxx
		${CXX} ${CXXFLAGS} ${INC} -MM $^ -MT $(patsubst %.d,%.o,$@) -MF $@ 

-include ${DEP}
