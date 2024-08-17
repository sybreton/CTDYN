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

$(EXE): dyna.o main.o profiles.o func_flow.o plg.o sort.o zbrent.o write_outputs.o bessel.o hunt.o \
	ccov.o cvect.o cio.o b2func.o b3func.o b4func.o b5func.o c1s2func.o c3s2func.o cdfunc.o
	$(LINKER) $(LDFLAGS) -o $@ $^  $(LIB) 

cio.mod: cio.f90 
	$(FC) $(FFLAGS) $(INCDIR) -c $< 
cvect.mod: cvect.f90 cio.mod 
	$(FC) $(FFLAGS) $(INCDIR) -c $< 
ccov.mod: ccov.f90 cio.mod
	$(FC) $(FFLAGS) $(INCDIR) -c $< 
b2func.mod: b2func.f90 cio.mod 
	$(FC) $(FFLAGS) $(INCDIR) -c $< 
b3func.mod: b3func.f90 cio.mod
	$(FC) $(FFLAGS) $(INCDIR) -c $< 
b4func.mod: b4func.f90 cio.mod
	$(FC) $(FFLAGS) $(INCDIR) -c $< 
b5func.mod: b5func.f90 cio.mod
	$(FC) $(FFLAGS) $(INCDIR) -c $< 
c1s2func.mod: c1s2func.f90 cio.mod
	$(FC) $(FFLAGS) $(INCDIR) -c $< 
c3s2func.mod: c3s2func.f90 cio.mod
	$(FC) $(FFLAGS) $(INCDIR) -c $< 
cdfunc.mod: cdfunc.f90 cio.mod
	$(FC) $(FFLAGS) $(INCDIR) -c $< 
bessel.mod: bessel.f90 cio.mod
	$(FC) $(FFLAGS) $(INCDIR) -c $< 
write_outputs.mod: write_outputs.f90 cio.mod ccov.mod cdfunc.mod \
                   bessel.mod hunt.mod func_flow.mod profiles.mod \
                   plg.mod 
	$(FC) $(FFLAGS) $(INCDIR) -c $< 
dyna.mod: dyna.f90 cint.f90 ccov.mod cvect.mod cio.mod \
	  b2func.mod b3func.mod b4func.mod b5func.mod  \
          c1s2func.mod c3s2func.mod cdfunc.mod bessel.mod \
          write_outputs.mod func_flow.mod profiles.mod \
          sort.mod 
	$(FC) $(FFLAGS) $(INCDIR) -c $< 
profiles.mod: profiles.f90 cio.mod 
	$(FC) $(FFLAGS) $(INCDIR) -c $< 
func_flow.mod: func_flow.f90 cio.mod  
	$(FC) $(FFLAGS) $(INCDIR) -c $< 
plg.mod: plg.f90  cio.mod
	$(FC) $(FFLAGS) $(INCDIR) -c $< 
zbrent.mod: zbrent.f90 cio.mod dyna.mod 
	$(FC) $(FFLAGS) $(INCDIR) -c $< 
hunt.mod: hunt.f90 cio.mod
	$(FC) $(FFLAGS) $(INCDIR) -c $< 
sort.mod: sort.f90 cio.mod
	$(FC) $(FFLAGS) $(INCDIR) -c $< 

cio.o: cio.f90 
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
cvect.o: cvect.f90 cio.mod 
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
ccov.o: ccov.f90 cio.mod
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
b2func.o: b2func.f90 cio.mod
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
b3func.o: b3func.f90 cio.mod
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
b4func.o: b4func.f90 cio.mod
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
b5func.o: b5func.f90 cio.mod
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
c1s2func.o: c1s2func.f90 cio.mod
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
c3s2func.o: c3s2func.f90 cio.mod
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
cdfunc.o: cdfunc.f90 cio.mod
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@ 
bessel.o: bessel.f90  cio.mod
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
write_outputs.o: write_outputs.f90 cio.mod ccov.mod cdfunc.mod \
	         bessel.mod 
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@

dyna.o: dyna.f90 cint.f90 ccov.mod cvect.mod cio.mod \
	b2func.mod b3func.mod b4func.mod b5func.mod  \
        c1s2func.mod c3s2func.mod cdfunc.mod bessel.mod \
        write_outputs.mod func_flow.mod profiles.mod \
        sort.mod
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
main.o: main.f90 cint.f90 cio.mod dyna.mod zbrent.mod 
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
profiles.o: profiles.f90 cio.mod 
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
func_flow.o: func_flow.f90 cio.mod  
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
plg.o: plg.f90  cio.mod
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
zbrent.o: zbrent.f90 cio.mod dyna.mod 
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
hunt.o: hunt.f90 cio.mod
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
sort.o: sort.f90 cio.mod
	$(FC) $(FFLAGS) $(INCDIR) -c $< -o $@
clean:	
clean:	
	rm *.o *.mod
	rm $(EXE)
