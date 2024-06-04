
FC = gfortran
FC90 = gfortran 
CPP = cpp
FPP = gfortran -E
LINKER = $(FC90)

#cobol
#FFLAGS =     -132 -O3 -vec_report0  -mcmodel=large -msse3
#muon
#FFLAGS =     -132 -O3  -mcmodel=large -msse3
#gollum (OSX)
#INTEL
#FFLAGS =      -132 -mkl -r8 -i4 -fast
#FFLAGS =      -132 -mkl -r8 -i4 -fast
#FFLAGS =      -132 -mkl -r8 -i4 -mcmodel=medium
#FFLAGS =      -132 -mkl -r8 -i4 -O3

#FFLAGS= -O3 -w -ffixed-line-length-132  -fdefault-integer-8 -fno-range-check -fdefault-double-8 -fdefault-real-8
FFLAGS= -O3 -ffixed-line-length-132  -fdefault-integer-8 -fno-range-check -fdefault-double-8 -fdefault-real-8

#FFLAGS =      -132 
#FFLAGS =    -132 -mp -w -vec_report0  -msss3  
#FFLAGS =  -O3  -Mextend -pc=80  -Mr8intrinsics    -Mcache_align 
#FFLAGS =  -O3  -Mr8 -Mextend -pc=80  -Mr8intrinsics  -Mcache_align -I/opt/mpich-pg64/include/
#F77FLAGS =   -r8 -i4 
#F77FLAGS =    
F77FLAGS =    
F90FLAGS = -I$(OBJDIR) -module $(OBJDIR) 
CFLAGS = 

#INTEL 
#LIB =  -lmkl_core -lmkl_intel_thread -lmkl_sequential -liomp5 -lpthread -lmkl_intel_lp64
#LIB =  -lmkl_core -lmkl_intel_thread -lmkl_sequential -liomp5  -lmkl_intel_lp64
#LIB =  -lmkl_core -lmkl_sequential -liomp5  -lmkl_intel_lp64
#LIB =  -lmkl_core -lmkl_sequential  -lmkl_intel_lp64
#LIB =  -lmkl_core -lmkl_intel_thread -lmkl_sequential -liomp5 -lmkl_intel_lp64
#LIB =  -lmkl_core -lmkl_sequential -liomp5 -lmkl_intel_lp64
#LIB =  -lmkl_core  -liomp5 -lmkl_intel_lp64
#LIB = -lmkl_core -lmkl_intel_thread -lmkl_sequential -liomp5 -lpthread -lmkl_intel_lp64
#LIB=-L/home/alfio/mesasdk/math-slots/crmath/lib  -lblas -llapack
#LIB=-L/home/alfio/lib/lapack-3.9.0 -llapack -lrefblas -lgfortran
#LIB=-L/home/alfio/lib/lapackfPIC -llapack -lrefblas -lgfortran


#@compilers to find lapack you may need to set:
#	   export LDFLAGS="-L/opt/homebrew/opt/lapack/lib"
#	        export CPPFLAGS="-I/opt/homebrew/opt/lapack/include"


#LIB =  -lmkl_core -lmkl_sequential  -lmkl_intel_lp64
#LIB =  -lmkl_core -lmkl_sequential  -lmkl_intel_lp64 -L/opt/intel/lib
#LIB =  -L/opt/intel/lib 
#LIB= -llapack64 -lblas64 -lgfortran
#LIB= -llapack -lblas  -L/opt/homebrew/opt/lapack/lib
#LIB= -llapack -lblas  -L/Applications/mesasdk/math-slots/crmath/lib
LIB= -llapack -lblas -lgfortran
LIB= -lblas -llapack -L/Users/alfio/lapack-3.12.0 -L/Users/alfio/BLAS-3.12.0



#cobol
#LDFLAGS = -i_dynamic
#muon
#DFLAGS = -Bdynamic
#DFLAGS = -i_dynamic
#osx
#LDFLAGS = -i_dynamic
CPPFLAGS = -DIF90_ -Dlapack
CPPFLAGS = -DIF90_ 
FCPPFL =  -fpp -DIF90_ -Dlapack
FCPPFL =  -fpp -DIF90_ 
OBJDIR = OBJECTS.Linux.PIV
INCDIR = -I.
PREDIR = PREDIR
EXE = ./spectra.x

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
