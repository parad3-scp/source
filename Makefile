CC=gcc
ICC=icc

LIBRARIES = -llapack -I/usr/include/lapacke/ -lm
LIBRARIES_INTEL = $(LIBRARIES) -L/opt/intel/compilers_and_libraries_2017.1.132/linux/compiler/lib/intel64_lin/
LIBRARIES_GCC = $(LIBRARIES)
CFLAGS=-I. -O3 -std=gnu99 -g -pg $(LIBRARIES_GCC)
CFLAGSP=-I. -O3 -std=gnu99 -g -pg $(LIBRARIES_GCC) -fopenmp
ICCFLAGS=-I.  -O3 -std=gnu99 -qopt-report5 -qopt-report-phase:openmp -qopenmp


ad3: ad3new.c ad3solve.c dynamicfd.c matrixinv.c util.c
	$(CC) -o ad3 ad3new.c ad3solve.c dynamicfd.c matrixinv.c util.c $(CFLAGS)


parad3: ad3new.c ad3solve.c dynamicfd.c matrixinv.c util.c
	$(CC) -o parad3 ad3new.c ad3solve.c dynamicfd.c matrixinv.c util.c $(CFLAGSP)

iad3: ad3new.c ad3solve.c dynamicfd.c matrixinv.c util.c
	$(ICC) -o iad3 ad3new.c ad3solve.c dynamicfd.c matrixinv.c util.c $(ICCFLAGS)

all: ad3 parad3 iad3

clean: 
	rm ad3 parad3 iad3
