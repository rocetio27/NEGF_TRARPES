
include ./make.inc

# modules
SRC_modules = modconstants.f90 modmain.f90 modpulse.f90 $(SRC_MPI) modmpi.f90 $(SRC_OMP) modomp.f90

# Elk program
SRC_main = main.f90

# main subroutines and functions
SRC_sub = \
 get_atom_pos.f90 genkline.f90 genhk.f90 \
 genafield_pump.f90 genafield_probe.f90 gengrid1d.f90 \
 nf.f90 phi_p.f90 hpp.f90 pw_coeff_bloch.f90\
 genhk_pw_basis_and_u_trans.f90\
 evolve_theta_arg_vec.f90 init_theta_arg_vec.f90

SRC = $(SRC_modules) $(SRC_main) $(SRC_sub)
OBJ = $(SRC:.f90=.o)
EXE = pes

#-------------------------------------------------------------------------------
# Suffix rules
#-------------------------------------------------------------------------------
.SUFFIXES: .o .f90
.f90.o:
	$(F90) $(F90_OPTS) -c $<

#-------------------------------------------------------------------------------
# Targets
#-------------------------------------------------------------------------------
elk:	$(OBJ)
	$(F90) $(F90_OPTS) -o $(EXE) $(OBJ) -ffree-line-length-512 -llapack -lblas

all:	elk

clean:
	rm -f *.o *.mod $(EXE)

