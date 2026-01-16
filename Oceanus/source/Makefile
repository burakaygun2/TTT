FC = gfortran-12

FFLAGS = -J ../mod \
    -I/usr/local/include -I. \
    -L/usr/local/lib \
    -fdefault-real-8 \
    -ffree-line-length-none \
    -g -fopenmp -fbacktrace -fcheck=all -fbounds-check -ffpe-trap=invalid,zero,overflow -Wall -Wextra -pedantic

FFLAGS_fast = -J ../mod \
    -I/usr/local/include -I. \
    -L/usr/local/lib \
    -fdefault-real-8 \
    -ffree-line-length-none \
    -Ofast -ffast-math -fopenmp -flto=auto -fwhole-program -march=native

SRC =./core/mod_parameters.f90
SRC+=./core/mod_variables.f90
SRC+=./core/mod_plib.f90
SRC+=./core/mod_allocation.f90
SRC+=./io/mod_read_f.f90
SRC+=./physics/mod_potential.f90
SRC+=./io/mod_write_f.f90
SRC+=./io/mod_write_h5.f90
SRC+=./core/mod_preprocess.f90
SRC+=./core/mod_preprocess_lsingle.f90
SRC+=./solver/mod_LU.f90
SRC+=./physics/mod_SHTns.f90
SRC+=./solver/mod_matrix_coeffs.f90
SRC+=./core/mod_math.f90
SRC+=./solver/mod_matrix.f90
# SRC+=./solver/mod_matrix_LPCK.f90
SRC+=./physics/mod_displacement.f90
SRC+=./solver/mod_rhs.f90
SRC+=./solver/mod_solver.f90
SRC+=./physics/mod_norms.f90
SRC+=./physics/mod_coriolis.f90
SRC+=./time/mod_time_integration.f90
SRC+=main.f90

LIBS = ../SHTns/libshtns.a -lfftw3 -lm -lc -llapack -lopenblas -I/usr/include/hdf5/serial \
	-L/usr/lib/x86_64-linux-gnu/hdf5/serial \
	-lhdf5 -lhdf5_fortran  \

TARGET = ../bin/oceanus.out

all: $(TARGET)

$(TARGET): $(SRC)
	$(FC) $(FFLAGS_fast) $(SRC) $(LIBS) -o $(TARGET)

clean:
	rm -f $(TARGET) ../mod/*.mod ./core/*.mod ./io/*.mod ./solver/*.mod

