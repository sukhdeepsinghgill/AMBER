target:	amber amber_v2 maps tidy clean

CC=gcc
CCM=mpicc
CFLAGS=-g
CFLAGSO=-O3 -ffast-math -march=native
ADDFLG=-flto -funroll-loops -fno-trapping-math -mavx2
########## change the following paths according to your setup ##########

FFTW_LIB=-L/path/to/local_install/lib/
GSL_LIB=-L/path/to/local_install/lib/
CFIO_LIB=-L/path/to/cfitsio/lib/

FFTW_INCLUDE=-I/path/to/local_install/include/
GSL_INCLUDE=-I/path/to/local_install/include/
CFIO_INCLUDE=-I/path/to/cfitsio/include/

#####################################################


amber:	amber.c vars.c amber_funcs.c
	$(CC) $(CFLAGSO) $(ADDFLG) $(FFTW_INCLUDE) $(CFIO_INCLUDE) $(GSL_INCLUDE) amber.c vars.c amber_funcs.c -fopenmp -lm $(FFTW_LIB) -lfftw3 -lfftw3_omp $(CFIO_LIB) -lcfitsio $(GSL_LIB) -lgsl -lgslcblas -o amber
	rm -rf *~ *.o	
	
	
mabs_final_v2:	amber_v2.c vars.c amber_funcs.c
	$(CC) $(CFLAGSO) $(ADDFLG) $(FFTW_INCLUDE) $(CFIO_INCLUDE) $(GSL_INCLUDE) amber_v2.c vars.c amber_funcs.c -fopenmp -lm $(FFTW_LIB) -lfftw3 -lfftw3_omp $(CFIO_LIB) -lcfitsio $(GSL_LIB) -lgsl -lgslcblas -o amber_v2
	rm -rf *~ *.o	


## The latest code is mabs_final and mabs_final_v2(with I). If confusion arises in result, match with mabs_fast_float (it is ultimate authority).
	
maps:	maps.c amber_funcs.c
	$(CC) $(CFLAGSO) $(ADDFLG) $(FFTW_INCLUDE) $(CFIO_INCLUDE) $(GSL_INCLUDE) maps.c mabs_funcs.c -fopenmp -lm $(FFTW_LIB) -lfftw3 -lfftw3_omp $(CFIO_LIB) -lcfitsio $(GSL_LIB) -lgsl -lgslcblas -o maps
	rm -rf *~ *.o
        
clean:
	rm -rf amber maps amber_v2 masked_* *~ *.o 
        
tidy:
	rm -rf masked_* modes *~ *.o bsout* apsout*
	
	
