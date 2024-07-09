FC = gfortran
FC90 = gfortran 
CPP = cpp
FPP = gfortran -E
LINKER = $(FC90)

FFLAGS = -O3 -ffixed-line-length-132  -fdefault-integer-8 -fno-range-check -fdefault-double-8 -fdefault-real-8
F77FLAGS =    
F90FLAGS = -I$(OBJDIR) -module $(OBJDIR) 
CFLAGS = 

LIB = -llapack -lblas -lgfortran
CPPFLAGS = -DIF90_ -Dlapack
CPPFLAGS = -DIF90_ 
FCPPFL =  -fpp -DIF90_ -Dlapack
FCPPFL =  -fpp -DIF90_ 
OBJDIR = OBJECTS.Linux.PIV
INCDIR = -I.
PREDIR = PREDIR
EXE = ../../ctdyn

vpath %.o .:$(OBJDIR)

$(EXE): dyna.o advection.o err.o func_flow.o plg.o  sort2.o zbrent.o writefield.o bessel.o hunt.o
	$(LINKER) $(LDFLAGS) -o $@ $^  $(LIB) 

dyna.o: dyna.f90 cint.f90 ccov.f90 cvect.f90 cio.f90
	$(FC) $(FFLAGS) $(INCDIR) $(F77FLAGS) -c $< -o $@
advection.o: advection.f90 cint.f90  
	$(FC) $(FFLAGS) $(INCDIR) $(F77FLAGS) -c $< -o $@
err.o: err.f90 cio.f90 
	$(FC) $(FFLAGS) $(INCDIR) $(F77FLAGS) -c $< -o $@
func_flow.o: func_flow.f90 b2func.f90 b3func.f90 b4func.f90 b5func.f90 c1s2func.f90 c3s2func.f90 cdfunc.f90 cio.f90  
	$(FC) $(FFLAGS) $(INCDIR) $(F77FLAGS) -c $< -o $@
plg.o: plg.f90  
	$(FC) $(FFLAGS) $(INCDIR) $(F77FLAGS) -c $< -o $@
zbrent.o: zbrent.f90  
	$(FC) $(FFLAGS) $(INCDIR) $(F77FLAGS) -c $< -o $@
outw.o: writefield.f90 cio.f90 ccov.f90 load_tor.f90 load_pol.f90 write_poloidal.f90 write_toroidal.f90 load_vec.f90
	$(FC) $(FFLAGS) $(INCDIR) $(F77FLAGS) -c $< -o $@
bessel.o: bessel.f90 
	$(FC) $(FFLAGS) $(INCDIR) $(F77FLAGS) -c $< -o $@
hunt.o: hunt.f90 
	$(FC) $(FFLAGS) $(INCDIR) $(F77FLAGS) -c $< -o $@
clean:	
clean:	
	rm *.o
	rm $(EXE)
