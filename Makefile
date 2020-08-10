#=========================== Makefile for Poseidon ============================#
#
#

#NEWTON_DIR = ./Code/Newton
#NEWTON_DRV = ./Code/Newton/Drivers
#NEWTON_SRC = ./Code/Newton/Code
#NEWTON_OBJ = ./Code/Newton/Obj
#NEWTON_BIN = ./Code/Newton/Bin

COM_SRC = ./Code/COMMON
DRV_SRC = ./Code/Newton/Drivers
SRC	= ./Code/Newton
OBJ	= ./Obj
BIN	= ./Bin



COMP	= gfortran
#FLAGs	= -g
#FLAGs   =  -O2
#FLAGs   =  -g -fcheck=all
FLAGs   =  -g -fbacktrace -ffpe-trap=invalid,zero -ffpe-summary=invalid,zero -fbounds-check
LFLAGs	= -framework Accelerate
EXT	= f90

#------------------------------ Set Source Code -------------------------------#



PROG	= Poseidon
INPUT	= 
OUTPUT	=



#--------------------------  List Code Componenets ----------------------------#



CODE_com = com.constants.o			 \
	   com.io.o


CODE_o  = z.Global_Variables_And_Parameters.o	 \
	  z.IO_Functions_Module.o		 \
	  z.Additional_Functions_Module.o	 \
	  z.Cholesky_Module.o			 \
	  z.Mesh_Module.o			 \
	  z.Stiffness_Matrix_Module.o		 \
	  z.Source_Vector_Module.o    		 \
	  z.Boundary_Conditions_Module.o  	 \
	  z.Test_Functions_Module.o		 \
	  z.Linear_Solvers_And_Preconditioners.o \
	  z.Coefficient_Vector_Module.o		 \
	  z.Poseidon_Main_Module.o		 \
	  z.Error_Analysis_Functions.o		 



#---------------------------- Compile Code ------------------------------------#

$(CODE_com):%.o: $(COM_SRC)/%.$(EXT)
	$(COMP) -J$(OBJ) -I$(OBJ) $(FLAGs) -c $< -o $(OBJ)/$@

$(CODE_o):%.o: $(SRC)/%.$(EXT)
	$(COMP) -J$(OBJ) -I$(OBJ) $(FLAGs) -c $< -o $(OBJ)/$@







# Core Program 


Poseidon: $(CODE_com) $(CODE_o)


PoseidonLib: $(CODE_com) $(CODE_o)
	ar crv $(OBJ)/poseidon.a $(OBJ)/*.o






# Program With Driver


main: $(CODE_com) $(CODE_o)
	@echo "         compiling with $(COMP) :" 
	$(COMP) -J$(OBJ) -I$(OBJ) $(FLAGs) -c $(DRV_SRC)/d.main.f90 -o $(OBJ)/d.main.o
	$(COMP) $(FLAGs) $(OBJ)/*.o  -o $(BIN)/main.x $(LFLAGs)
	@echo ">>> compiled on `hostname -s` with $(COMP) <<<"
	rm -f $(OBJ)/com.*.o
	rm -f $(OBJ)/d.*.o
	rm -f $(OBJ)/z.*.o




# Cholesky Timing Code


Chol_Time: $(CODE_com) $(CODE_o)
	@echo "         compiling with $(COMP) :"
	$(COMP) -J$(OBJ) -I$(OBJ) $(FLAGs) -c $(DRV_SRC)/d.Chol_Time.f90 -o $(OBJ)/d.Chol_Time.o
	$(COMP) $(FLAGs) $(OBJ)/*.o  -o $(BIN)/Chol_Time.x $(LFLAGs)
	@echo ">>> compiled on `hostname -s` with $(COMP) <<<"
	rm -f $(OBJ)/com.*.o
	rm -f $(OBJ)/d.*.o
	rm -f $(OBJ)/z.*.o






# Prints Inf and One error norms and convergence rate

Test_Converge: $(CODE_com) $(CODE_o)
	@echo "         compiling with $(COMP) :"
	$(COMP) -J$(OBJ) -I$(OBJ) $(FLAGs) -c $(DRV_SRC)/d.Test_Converge.f90 -o $(OBJ)/d.Test_Converge.o
	$(COMP) $(FLAGs) $(OBJ)/*.o  -o $(BIN)/Test_Converge.x $(LFLAGs)
	@echo ">>> compiled on `hostname -s` with $(COMP) <<<"
	rm -f $(OBJ)/com.*.o
	rm -f $(OBJ)/d.*.o
	rm -f $(OBJ)/z.*.o







# Prints Inf and One error norms and convergence rate

Test_Compare: $(CODE_com) $(CODE_o)
	@echo "         compiling with $(COMP) :"
	$(COMP) -J$(OBJ) -I$(OBJ) $(FLAGs) -c $(DRV_SRC)/d.Test_Compare.f90 -o $(OBJ)/d.Test_Compare.o
	$(COMP) $(FLAGs) $(OBJ)/*.o  -o $(BIN)/Test_Compare.x $(LFLAGs)
	@echo ">>> compiled on `hostname -s` with $(COMP) <<<"
	rm -f $(OBJ)/com.*.o
	rm -f $(OBJ)/z.*.o
	rm -f $(OBJ)/d.*.o







Test_Compare_Map: $(CODE_com) $(CODE_o)
	@echo "         compiling with $(COMP) :"
	$(COMP) -J$(OBJ) -I$(OBJ) $(FLAGs) -c $(DRV_SRC)/d.Test_Compare_Map.f90 -o $(OBJ)/d.Test_Compare_Map.o
	$(COMP) $(FLAGs) $(OBJ)/*.o  -o $(BIN)/Test_Compare_Map.x $(LFLAGs)
	@echo ">>> compiled on `hostname -s` with $(COMP) <<<"
	rm -f $(OBJ)/com.*.o
	rm -f $(OBJ)/z.*.o
	rm -f $(OBJ)/d.*.o






Test_L2_Error: $(CODE_com) $(CODE_o)
	@echo "         compiling with $(COMP) :"
	$(COMP) -J$(OBJ) -I$(OBJ) $(FLAGs) -c $(DRV_SRC)/d.Test_L2_Error.f90 -o $(OBJ)/d.Test_L2_Error.o
	$(COMP) $(FLAGs) $(OBJ)/*.o  -o $(BIN)/Test_L2_Error.x $(LFLAGs)
	@echo ">>> compiled on `hostname -s` with $(COMP) <<<"
	rm -f $(OBJ)/com.*.o
	rm -f $(OBJ)/z.*.o
	rm -f $(OBJ)/d.*.o




Test_Time: $(CODE_com) $(CODE_o)
	@echo "         compiling with $(COMP) :"
	$(COMP) -J$(OBJ) -I$(OBJ) $(FLAGs) -c $(DRV_SRC)/d.Test_Time.f90 -o $(OBJ)/d.Test_Time.o
	$(COMP) $(FLAGs) $(OBJ)/*.o  -o $(BIN)/Test_Time.x $(LFLAGs)
	@echo ">>> compiled on `hostname -s` with $(COMP) <<<"
	rm -f $(OBJ)/com.*.o
	rm -f $(OBJ)/z.*.o
	rm -f $(OBJ)/d.*.o








# Tests Legendre-Gauss-Lobatto Quadrature by calculating known integratal of x^2 on [-1,1].

Test_LGL: $(CODE_o)
	@echo "         compiling with $(COMP) :"
	$(COMP) -J$(OBJ) -I$(OBJ) $(FLAGs) -c $(DRV_SRC)/d.Test_LGL.f90 -o $(OBJ)/d.Test_LGL.o
	$(COMP) $(FLAGs) $(OBJ)/*.o  -o $(BIN)/Test_LGL.x $(LFLAGs)
	@echo ">>> compiled on `hostname -s` with $(COMP) <<<"
	rm -f $(OBJ)/com.*.o
	rm -f $(OBJ)/z.*.o
	rm -f $(OBJ)/d.*.o






Print_Matrix: $(CODE_o)
	@echo "         compiling with $(COMP) :"
	$(COMP) -J$(OBJ) -I$(OBJ) $(FLAGs) -c $(DRV_SRC)/d.Print_Matrix.f90 -o $(OBJ)/d.Print_Matrix.o
	$(COMP) $(FLAGs) $(OBJ)/*.o  -o $(BIN)/Print_Matrix.x $(LFLAGs)
	@echo ">>> compiled on `hostname -s` with $(COMP) <<<"
	rm -f $(OBJ)/com.*.o
	rm -f $(OBJ)/z.*.o
	rm -f $(OBJ)/d.*.o



Matrix_Details: $(CODE_o)
	@echo "         compiling with $(COMP) :"
	$(COMP) -J$(OBJ) -I$(OBJ) $(FLAGs) -c $(DRV_SRC)/d.Matrix_Details.f90 -o $(OBJ)/d.Matrix_Details.o
	$(COMP) $(FLAGs) $(OBJ)/*.o  -o $(BIN)/Matrix_Details.x $(LFLAGs)
	@echo ">>> compiled on `hostname -s` with $(COMP) <<<"
	rm -f $(OBJ)/com.*.o
	rm -f $(OBJ)/z.*.o
	rm -f $(OBJ)/d.*.o



#---------------------------- Execute Code ------------------------------------#


run_main: 
	./$(BIN)/main.x

run_Chol_Time:
	./$(BIN)/Chol_Time.x

run_Test_Converge:
	./$(BIN)/Test_Converge.x


run_Test_Compare:
	./$(BIN)/Test_Compare.x


run_Test_Compare_Map:
	./$(BIN)/Test_Compare_Map.x


run_Test_L2_Error:
	./$(BIN)/Test_L2_Error.x


run_Test_Time:
	./$(BIN)/Test_Time.x

run_Test_LGL:
	./$(BIN)/Test_LGL.x

run_Print_Matrix:
	./$(BIN)/Print_Matrix.x

run_Matrix_Details:
	./$(BIN)/Matrix_Details.x


#------------------------------- Clean Up ------------------------------------#

clean:
	
	rm -f $(OBJ)/com.*.o
	rm -f $(OBJ)/z.*.o
	rm -f $(OBJ)/d.*.o
	rm -f $(OBJ)/*.mod DONE


