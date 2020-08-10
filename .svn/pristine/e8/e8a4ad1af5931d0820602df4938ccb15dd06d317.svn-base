# Variable Definition 

NEWTON_DIR = ./Code/Newton
NEWTON_DRV = ./Code/Newton/Drivers
NEWTON_SRC = ./Code/Newton/Code
NEWTON_OBJ = ./Code/Newton/Obj
NEWTON_BIN = ./Code/Newton/Bin

CFA_DIR = ./Code/CFA
CFA_DRV = ./Code/CFA/Drivers
CFA_SRC = ./Code/CFA/Code
CFA_OBJ = ./Code/CFA/Obj
CFA_BIN = ./Code/CFA/Bin


SUBDIRS = $(NEWTON_DIR) $(CFA_DIR)


Newton: 
	@$(MAKE) -f Makefile.Newton $(ARG)


CFA: $(C_CALL)
	@ make -f Makefile.CFA  $(ARG)



 
