# Makefile for building a shared library (dylib) from Fortran and C code for VELUGA
# 	--- For LINUX
#

# Conda dependence
CONDA_INC = $(CONDA_PREFIX)/include
CONDA_LIB = $(CONDA_PREFIX)/lib

# Source Files
BASES	= find_domain jsrd_part_totnum jsrd_part get_ptcl prop_time jsamr2cell_totnum jsamr2cell rv_match prop_sfr js_kde_2dhisto js_kde_gauss_binning js_kde_gauss_pts js_gasmap get_hash get_merit get_merit2 js_indmatch get_contam js_getpt_ft

F90_SRCS	= $(BASES:=.f90)
F90_OBJS	= $(BASES:=.o)
C_SRCS		= $(BASES:=_c2f.c)
C_OBJS		= $(BASES:=_c2f.o)
DYLIBS		= $(BASES:=.so)

# Compilers
FC = gfortran
CC = gcc

# Flags
PIC  = -fPIC
OMP  = -fopenmp
LMEM = -mcmodel=large
NORANGE = -fno-range-check

# Compiles
all: ckdtree $(DYLIBS)

ckdtree:
	$(FC) $(OMP) $(PIC) -I$(CONDA_INC) -c js_kdtree.f90 -o js_kdtree.o

%.o: %.f90
	$(FC) $(OMP) $(PIC) $(NORANGE) -I$(CONDA_INC) -c $< -o $@

%_c2f.o: %_c2f.c
	$(CC) $(OMP) $(PIC) -I$(CONDA_INC) -c $< -o $@

%.so: %.o %_c2f.o
	$(FC) -shared $(OMP) -L$(CONDA_LIB) $^ js_kdtree.o -o $@


clean:
	rm -f *.o

.PHONT: all clean
