FC = gfortran

FFLAGS = -Ofast -fopenmp -funroll-loops -finline-functions \
         -fno-unsafe-math-optimizations -ffloat-store -march=native

SRC = types.f90 utils.f90 math.f90 materials.f90 components.f90 sources.f90 \
      base_inductors.f90 circular_inductors.f90 \
      inductors_kernel.f90

OBJ = $(SRC:.f90=.o)

MOD = types.mod utils_library.mod math_library.mod materials_library.mod sources_library.mod \
      components_library.mod base_inductors_library.mod \
      circular_inductors_library.mod \
      inductor_solvers.mod

# =========================
# OBJETIVO PRINCIPAL
# =========================

all: $(OBJ)

# =========================
# COMPILACIÓN DE MÓDULOS
# =========================

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

# =========================
# LIMPIEZA
# =========================

clean:
	rm -f $(OBJ) $(MOD)

