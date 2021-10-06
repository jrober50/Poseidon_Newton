POSEIDON_NEWTON_DIR ?= /ccsopen/home/jrober50/Poseidon_Newton

## Select Number of Processes (Mostly for MacBook)

NPROCS=1

## Select Machine

MACHINE_NAME    =Ascent


## Select Mode

CMODE   =DEBUG
#CMODE  =OPTIMIZE


## Compile with Openmp

#OPENMP_MODE    =ON
OPENMP_MODE     =OFF


#OPENACC_MODE    =ON
OPENACC_MODE    =OFF


## Compile with PETSc

#PETSC_MODE      =ON
PETSC_MODE     =OFF


## Compile with HDF5

HDF5_MODE       =ON
#HDF5_MODE      =OFF


## Compile with MPI

MPI_MODE        =ON
#MPI_MODE       =OFF

AMREX_MODE      =OFF


include $(POSEIDON_NEWTON_DIR)/Build/Makefile_Core

Main : $(POSEIDON_o)
	@echo "         compiling with $(COMP_$(MACHINE_NAME)) :"
	$(FORT)$(POSEIDON_DEFS) $(FLAGS) $(INCLUDE) -c Drivers/Main.f90  -o $(OBJ)/Driver.o $(OUTPUT)
	$(FORT) $(FLAGS) $(OBJ)/*.o -o $(BIN)/Poseidon_Newton.x
	@echo ">>> compiled on `hostname -s` with $(FORT_$(MACHINE_NAME)) <<<"


run : 
	./$(BIN)/Poseidon_Newton.x


clean:
	rm -f $(OBJ)/*.o
	rm -f $(OBJ)/*.mod
	rm -f $(BIN)/*.x 

clean_output:
	rm -f $(OUT)/*.out

