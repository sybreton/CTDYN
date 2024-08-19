FC = gfortran
FFLAGS = -O3 -ffixed-line-length-132 -fdefault-integer-8 -fno-range-check -fdefault-double-8 -fdefault-real-8
LIB = -llapack -lblas -lgfortran

SRCDIR = src/CTDYN/
BINDIR = bin/
OBJDIR = obj/
TARGET = ctdyn

# Retrieve list of the source files
SRC = $(wildcard $(addsuffix *.f90,$(SRCDIR)))
# Generate list of the object files
OBJ = $(addprefix $(OBJDIR), $(patsubst %.f90, %.o, $(notdir $(SRC))))
# Set directory where module will be written
MODDIR = $(OBJDIR)
FFLAGS += -J $(MODDIR)

VPATH = $(SRCDIR)

$(TARGET) : $(OBJ)
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
