FC = gfortran
FFLAGS = -O3 -ffixed-line-length-132 -fdefault-integer-8 -fno-range-check -fdefault-double-8 -fdefault-real-8
LIB = -llapack -lblas -lgfortran

SRCDIR = src/CTDYN/
BINDIR = bin/
OBJDIR = obj/
TARGET = ctdyn

# List of the source files
SRC_BASE = kind_parameters 
SRC_BASE += util cio 
SRC_BASE += b2func b3func b4func b5func c1s2func c3s2func plg 
SRC_BASE += bessel cdfunc stellar_profiles 
SRC_BASE += func_flow write_outputs dyna zbrent main
SRC = $(foreach item,$(SRC_BASE), $(SRCDIR)$(item).f90)

# Generate list of the object files
OBJ = $(addprefix $(OBJDIR), $(patsubst %.f90, %.o, $(notdir $(SRC))))
# Set directory where module will be written
MODDIR = $(OBJDIR)
FFLAGS += -J $(MODDIR)

VPATH = $(SRCDIR)

$(TARGET) : $(OBJ)
	@echo $(SRC_BASE)
	@echo $(SRC)
	@echo Linking...
	@mkdir -p $(BINDIR)
	@$(FC) $(FFLAGS) -o $(BINDIR)$@ $(OBJ) $(LIB)

$(OBJDIR)%.o : %.f90
	@echo Compiling $< in $@...
	@mkdir -p $(OBJDIR)
	@$(FC) $(FFLAGS) -c -o $@ $^

clean :
	@$(RM) -r $(OBJDIR)
	@$(RM) -r $(BINDIR)
