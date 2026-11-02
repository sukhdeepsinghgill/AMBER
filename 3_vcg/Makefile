CC=gcc
CFLAGS=-c
CFLAGSO=-O3 -ffast-math -march=native
ADDFLG=-flto -funroll-loops -mavx2

FFTW_LIB=-L/scratch/sukhdeep/project_MABS/libraries/local_install/lib/
GSL_LIB=-L/scratch/sukhdeep/project_MABS/libraries/local_install/lib/
CFITS_LIB=-L/scratch/sukhdeep/project_MABS/libraries/cfitsio/lib/

FFTW_INCLUDE=-I/scratch/sukhdeep/project_MABS/libraries/local_install/include/
GSL_INCLUDE=-I/scratch/sukhdeep/project_MABS/libraries/local_install/include/
CFITS_INCLUDE=-I/scratch/sukhdeep/project_MABS/libraries/cfitsio/include/


	
read_fits_funcmfs.o:	read_fits_funcmfs.c
	$(CC) $(CFLAGS) $(FFTW_INCLUDE) $(GSL_INCLUDE) $(CFITS_INCLUDE) $(CFITS_LIB) $(FFTW_LIB) read_fits_funcmfs.c -o read_fits_funcmfs.o

tgefuncs.o:	tgefuncs.c
	$(CC) $(CFLAGS) $(FFTW_INCLUDE) $(GSL_INCLUDE) $(CFITS_INCLUDE) $(CFITS_LIB) $(FFTW_LIB) tgefuncs.c -o tgefuncs.o
		

main.o:	main.c
	$(CC) $(CFLAGS) $(FFTW_INCLUDE) $(GSL_INCLUDE) $(CFITS_INCLUDE) $(CFITS_LIB) $(FFTW_LIB) main.c -o main.o

vcg:		read_fits_funcmfs.o tgefuncs.o main.o
	$(CC) $(CFLAGSO) $(ADDFLG) read_fits_funcmfs.o tgefuncs.o main.o $(FFTW_INCLUDE) $(GSL_INCLUDE) $(CFITS_INCLUDE) $(CFITS_LIB) -lm -lcfitsio $(FFTW_LIB) -lfftw3 -lgsl -lgslcblas -o vcg

tidy:
	rm -rf *.o
	rm -rf *~
	
clean:
	rm -rf vcg tge
	rm -rf *.o
	rm -rf *~
	


