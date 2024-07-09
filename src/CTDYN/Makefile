FC = gfortran
FC90 = gfortran 
CPP = cpp
FPP = gfortran -E
LINKER = $(FC90)

FFLAGS = -O3 -ffixed-line-length-132  -fdefault-integer-8 -fno-range-check -fdefault-double-8 -fdefault-real-8
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

$(EXE): dyna.o advection.o err.o func_flow.o plg.o sort2.o zbrent.o writefield.o bessel.o hunt.o \
	ccov.o cvect.o cio.o b2func.o b3func.o b4func.o b5func.o c1s2func.o c3s2func.o cdfunc.o
	$(LINKER) $(LDFLAGS) -o $@ $^  $(LIB) 

cio.mod: cio.f90 
	$(FC) $(FFLAGS) $(INCDIR) -c $< 
cvect.mod: cvect.f90 cio.mod 
	$(FC) $(FFLAGS) $(INCDIR) -c $< 
ccov.mod: ccov.f90 cio.mod
	$(FC) $(FFLAGS) $(INCDIR) -c $< 
b2func.mod: b2func.f90 
	$(FC) $(FFLAGS) $(INCDIR) -c $< 
b3func.mod: b3func.f90 
	$(FC) $(FFLAGS) $(INCDIR) -c $< 
b4func.mod: b4func.f90 
	$(FC) $(FFLAGS) $(INCDIR) -c $< 
b5func.mod: b5func.f90 
	$(FC) $(FFLAGS) $(INCDIR) -c $< 
c1s2func.mod: c1s2func.f90 
	$(FC) $(FFLAGS) $(INCDIR) -c $< 
c3s2func.mod: c3s2func.f90 
	$(FC) $(FFLAGS) $(INCDIR) -c $< 
cdfunc.mod: cdfunc.f90 
	$(FC) $(FFLAGS) $(INCDIR) -c $< 
bessel.mod: bessel.f90 
	$(FC) $(FFLAGS) $(INCDIR) -c $< 

cio.o: cio.f90 
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
cvect.o: cvect.f90 cio.mod 
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
ccov.o: ccov.f90 cio.mod
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
b2func.o: b2func.f90 
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
b3func.o: b3func.f90 
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
b4func.o: b4func.f90 
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
b5func.o: b5func.f90 
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
c1s2func.o: c1s2func.f90 
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
c3s2func.o: c3s2func.f90 
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
cdfunc.o: cdfunc.f90 
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@ 
bessel.o: bessel.f90  
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@

dyna.o: dyna.f90 cint.f90 ccov.mod cvect.mod cio.mod \
	b2func.mod b3func.mod b4func.mod b5func.mod  \
        c1s2func.mod c3s2func.mod cdfunc.mod bessel.mod 
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
advection.o: advection.f90 cint.f90  
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
err.o: err.f90 cio.mod 
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
func_flow.o: func_flow.f90 cio.mod  
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
plg.o: plg.f90  
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
zbrent.o: zbrent.f90  
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
writefield.o: writefield.f90 cio.mod ccov.mod cdfunc.mod \
	bessel.mod load_tor.f90 load_pol.f90 write_poloidal.f90 \
	write_toroidal.f90 load_vec.f90
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
hunt.o: hunt.f90 
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
sort2.o: sort2.f90 
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
clean:	
clean:	
	rm *.o *.mod
	rm $(EXE)
