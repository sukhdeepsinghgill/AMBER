CC = gcc
CFLAGS ?= -O2 -std=gnu11 -Wall -Wextra -fopenmp
LDLIBS = -lfftw3_omp -lfftw3 -lgsl -lgslcblas -lcfitsio -lm

.PHONY: all clean
all: 3_vcg/vcg 4_fbe/amber 4_fbe/amber_v2

3_vcg/vcg: 3_vcg/read_fits_funcmfs.c 3_vcg/tgefuncs.c 3_vcg/main.c 3_vcg/tgefuncs.h 3_vcg/read_fits_func.h
	$(CC) $(CPPFLAGS) $(CFLAGS) 3_vcg/read_fits_funcmfs.c 3_vcg/tgefuncs.c 3_vcg/main.c $(LDFLAGS) $(LDLIBS) -o $@

4_fbe/amber: 4_fbe/amber.c 4_fbe/vars.c 4_fbe/amber_funcs.c 4_fbe/vars.h 4_fbe/mabs_funcs.h
	$(CC) $(CPPFLAGS) $(CFLAGS) 4_fbe/amber.c 4_fbe/vars.c 4_fbe/amber_funcs.c $(LDFLAGS) $(LDLIBS) -o $@

4_fbe/amber_v2: 4_fbe/amber_v2.c 4_fbe/vars.c 4_fbe/amber_funcs.c 4_fbe/vars.h 4_fbe/mabs_funcs.h
	$(CC) $(CPPFLAGS) $(CFLAGS) 4_fbe/amber_v2.c 4_fbe/vars.c 4_fbe/amber_funcs.c $(LDFLAGS) $(LDLIBS) -o $@

clean:
	rm -f 3_vcg/vcg 4_fbe/amber 4_fbe/amber_v2
