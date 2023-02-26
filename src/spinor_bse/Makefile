SHELL = /bin/sh

MKLROOT = /opt/intel/compilers_and_libraries/linux/mkl
INTELLINK = ${MKLROOT}/lib/intel64/libmkl_scalapack_ilp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_ilp64.a -Wl,--end-group -liomp5 -lpthread -lm -ldl

#Linux_LIB = -shared-intel -Wl,--start-group /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_sequential.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_core.a -lfftw3 -fopenmp
#Linux_LIB = -shared-intel -Wl,--start-group /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_intel_ilp64.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_lapack95_ilp64.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_blas95_ilp64.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_intel_thread.a /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64/libmkl_core.a -qopenmp -lpthread -lfftw3 

Linux_LIB =  -lfftw3 ${INTELLINK}
Linux_MYLIB = 
Linux_CLIB = -lm 
Linux_FCLIB = 	

CYGWIN_NT-6.3_LIB = -llapack -lblas -lfftw3_threads -lfftw3 -lpthread
CYGWIN_NT-6.3_MYLIB = 
CYGWIN_NT-6.3_CLIB =  -lm
CYGWIN_NT-6.3_FCLIB = 	

Windows_NT_LIB = -llapack -lblas -lfftw3_threads -lfftw3 -lpthread
Windows_NT_MYLIB = 
Windows_NT_CLIB =  -lm
Windows_NT_FCLIB = 	

LIB = $(${OS}_LIB) $(${OS}_CLIB) $(${OS}_FCLIB) $(${OS}_MYLIB)

# flags ...
Linux_OPTFLAGS = -DMKL_ILP64 -O3 -xSSE3 -fp-model=fast -qopenmp -static 
#Linux_OPTFLAGS = -O0 -g -qopenmp
Linux_CFLAGS = $(${OS}_OPTFLAGS) -DFFT_FFTW
Linux_FFLAGS = $(${OS}_OPTFLAGS)

CYGWIN_NT-6.3_OPTFLAGS = -O3 -ffast-math
CYGWIN_NT-6.3_CFLAGS = $(${OS}_OPTFLAGS) -DFFT_FFTW
CYGWIN_NT-6.3_FFLAGS = $(${OS}_OPTFLAGS) 

Windows_NT_OPTFLAGS = -O3 -ffast-math
Windows_NT_CFLAGS = $(${OS}_OPTFLAGS) -DFFT_FFTW
Windows_NT_FFLAGS = $(${OS}_OPTFLAGS) 


MAINNAM = bs

#compiler ...
Linux_CC = icx
Linux_FF = ifx
Linux_LD = icx

CYGWIN_NT-6.3_CC = gcc
CYGWIN_NT-6.3_FF = gfortran
CYGWIN_NT-6.3_LD = gcc

Windows_NT_CC = gcc
Windows_NT_FF = gfortran
Windows_NT_LD = gcc

OBJECTS = \
	main.o init.o size.o norm.o nerror.o read.o hartree.o single.o energy.o interpolate.o \
	bethe-salpeter.o diag.o rand.o dipole.o hamiltonian.o write.o angular.o



# compilation ...

.f.o:
	$(${OS}_FF) $(${OS}_FFLAGS) -c  $*.f
.c.o:
	$(${OS}_CC) $(${OS}_CFLAGS) -c  $*.c  

$(MAINNAM): $(OBJECTS) 
	$(${OS}_LD) -o $(MAINNAM).x $(${OS}_CFLAGS) $(OBJECTS) $(LIB)

convert: convert.o
	$(${OS}_CC) $(${OS}_CFLAGS) -o convert.x convert.o

clean:
	/bin/rm *.o *.x
