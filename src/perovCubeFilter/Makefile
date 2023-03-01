SHELL = /bin/sh

#MKLROOT = /opt/intel/compilers_and_libraries/linux/mkl
INTELLINK =  ${MKLROOT}/lib/intel64/libmkl_scalapack_ilp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_ilp64.a -Wl,--end-group -liomp5 -lpthread -lm -ldl

Linux_LIB =  ${INTELLINK} -L${HOME}/usr/lib -lfftw3
Linux_MYLIB = 
Linux_CLIB = -lm
Linux_FCLIB = 	

CYGWIN_NT-6.3_LIB = -llapack -lblas -lfftw3_threads -lfftw3 -lpthread
CYGWIN_NT-6.3_MYLIB = 
CYGWIN_NT-6.3_CLIB =  -lm
CYGWIN_NT-6.3_FCLIB = 	

LIB = $(${OS}_LIB) $(${OS}_CLIB) $(${OS}_FCLIB) $(${OS}_MYLIB)

# flags ...
#Linux_OPTFLAGS = -fast -qopenmp -static
Linux_OPTFLAGS = -DMKL_ILP64 -O3 -xSSE3  -qopenmp -static
#Linux_OPTFLAGS = -O0 -g
Linux_CFLAGS = $(${OS}_OPTFLAGS) -DFFT_FFTW -I${HOME}/usr/include -I /opt/rocks/include -I"${MKLROOT}/include" -Wall -pedantic 
Linux_FFLAGS = $(${OS}_OPTFLAGS)

CYGWIN_NT-6.3_OPTFLAGS = -O3 -ffast-math 
CYGWIN_NT-6.3_CFLAGS = $(${OS}_OPTFLAGS) -DFFT_FFTW
CYGWIN_NT-6.3_FFLAGS = $(${OS}_OPTFLAGS) 

MAINNAM = spinOrbitFilter

#compiler ...
Linux_CC = icx
Linux_FF = ifx
Linux_LD = icx
CYGWIN_NT-6.3_CC = gcc
CYGWIN_NT-6.3_FF = gfortran
CYGWIN_NT-6.3_LD = gcc

OBJECTS = \
	main.o init.o size.o norm.o nerror.o coeff.o read.o interpolate.o \
	ortho.o energy.o hamiltonian.o hnorm.o filter.o Hmat.o dipole.o rand.o \
	projectors.o angular.o write.o

CUBEOBJECTS = \
		makecube.o write.o size.o init.o nerror.o read.o rand.o interpolate.o projectors.o norm.o
# compilation ...

.f.o:
	$(${OS}_FF) $(${OS}_FFLAGS) -c  $*.f
.c.o:
	$(${OS}_CC) -DOS_$(OS) $(${OS}_CFLAGS) -c  $*.c  


$(MAINNAM): $(OBJECTS) 
	$(${OS}_LD) -o $(MAINNAM).x $(${OS}_CFLAGS) $(OBJECTS) $(LIB)

all: $(MAINNAM) cube

cube: $(CUBEOBJECTS)
	$(${OS}_LD) -o makecube.x $(${OS}_CFLAGS) $(CUBEOBJECTS) $(LIB)

clean:
	/bin/rm *.o *.x
