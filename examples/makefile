# Compilador Fortran
FC = gfortran

# Flags de compilación
# -fopenmp: OpenMP
# -I../src: código fuente propio
# -I/usr/include: ruta donde está hdf5.mod
FFLAGS = -fopenmp -I../src -I/usr/include

# Flags de enlace
# Librerías HDF5 Fortran y C
LDFLAGS = -L/usr/lib -lhdf5_fortran -lhdf5

# Flags de Intel MKL
MKL_FLAGS = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

# Objetos a compilar
OBJ = ../src/*.o simulation.o

# Target principal
simulation: $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) -o simulation $(LDFLAGS) $(MKL_FLAGS)

# Compilar el archivo principal
simulation.o: simulation.f90
	$(FC) $(FFLAGS) -c simulation.f90

# Limpiar archivos objeto y ejecutable
clean:
	rm -f simulation.o simulation

