CC       = gcc
MPICC    = mpicc.mpich

OPENMP   = -fopenmp
CFLAGS   = -Wall -Wimplicit-function-declaration -O3 -std=c99 -march=native 

TARGETS = pi_mpi_collectives pi_seq pi_mpi pi_mpi_omp

all: $(TARGETS)

pi_mpi_collectives: pi_mpi.c
	$(MPICC) $< $(CFLAGS) -D_COLLECTIVES_ $(CINCL) -o $@ $(CLIBS)

pi_mpi: pi_mpi.c
	$(MPICC) $< $(CFLAGS) $(CINCL) -o $@ $(CLIBS)

pi_mpi_omp: pi_mpi_omp.c
	$(MPICC) $< $(CFLAGS) $(OPENMP) $(CINCL) -o $@ $(CLIBS)

pi_seq: pi_seq.c
	$(CC) $< $(CFLAGS) -o $@

clean:
	rm -rf $(TARGETS) *.o *.pe* *.po* *.o* *.e*
ultraclean:
	rm -rf $(TARGETS) *.o *.pe* *.po* *.o* *.e* *.txt *.ps
