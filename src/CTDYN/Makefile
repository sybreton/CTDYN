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
EXE = ../../spectra.x

vpath %.o .:$(OBJDIR)

$(EXE): dyna.o advection.o err.o func_flow.o plg.o  sort2.o zbrent.o writefield.o bessel.o hunt.o
	$(LINKER) $(LDFLAGS) -o $@ $^  $(LIB) 

dyna.o: dyna.f cint ccov cvect cio
	$(FC) $(FFLAGS) $(INCDIR) $(F77FLAGS) -c $< -o $@
advection.o: advection.f cint  
	$(FC) $(FFLAGS) $(INCDIR) $(F77FLAGS) -c $< -o $@
err.o: err.f cio 
	$(FC) $(FFLAGS) $(INCDIR) $(F77FLAGS) -c $< -o $@
func_flow.o: func_flow.f b2func b3func b4func b5func c1s2func c3s2func cdfunc cio  
	$(FC) $(FFLAGS) $(INCDIR) $(F77FLAGS) -c $< -o $@
plg.o: plg.f  
	$(FC) $(FFLAGS) $(INCDIR) $(F77FLAGS) -c $< -o $@
zbrent.o: zbrent.f  
	$(FC) $(FFLAGS) $(INCDIR) $(F77FLAGS) -c $< -o $@
outw.o: wwritefield.f cio ccov load_tor.f load_pol.f write_poloidal.f write_toroidal.f load_vec.f
	$(FC) $(FFLAGS) $(INCDIR) $(F77FLAGS) -c $< -o $@
bessel.o: bessel.f 
	$(FC) $(FFLAGS) $(INCDIR) $(F77FLAGS) -c $< -o $@
hunt.o: hunt.f 
	$(FC) $(FFLAGS) $(INCDIR) $(F77FLAGS) -c $< -o $@
clean:	
clean:	
	rm *.o
	rm $(EXE)
