FORMAT=ADD0US
FFTW_PATH = /home/p/pen/malvarez/fftw-2.1.5

F90         = mpif90
F77         = mpif77
CC          = mpicc
CC_COMPILE  = $(CC) -c
F90_OPT     = -O3
F90_COMPILE = $(F90_OPT) -c
F77_COMPILE = $(F90_OPT) -c


FFTLIB =  -L$(FFTW_PATH)/lib -lsrfftw -lsfftw -lsrfftw_mpi -lsfftw_mpi

.SUFFIXES: .f .f90 .F .F90

%.o : $.mod
	$(F90) $(F90_COMPILE) $< 

%.o: %.f90 
	$(F90) $(F90_COMPILE) -I$(FFTW_PATH)/include  $< 

%.o: %.f 
	$(F77) $(F77_COMPILE) -I$(FFTW_PATH)/include  $< 

%.o: %.F
	$(F77) $(F90_COMPILE) -D$(FORMAT) -I$(FFTW_PATH)/include $< 

%.o: %.F90
	$(F90) $(F90_COMPILE) -D$(FORMAT) -I$(FFTW_PATH)/include $< 

exec = powerspectra
objs = grafic_types.o grid.o  \
	   transform.o textlib.o timing_diagnostics.o \
	   mpivars.o memory_management.o mpnorm.o \
	   pk_fftw.o correlation.o pktable.o pks2grid.o powerspectra.o 

all_execs = $(exec)


$(exec): $(objs)
	$(F90) $(F90_OPT) $(objs) $(FFTLIB) -o $(exec)
	cp $(exec) $(HOME)/bin

clean:
	rm -f *.o *.mod $(all_execs)

