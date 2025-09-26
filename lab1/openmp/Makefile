CC = gcc
OPENMP	 = 
CFLAGS   = -Wall -O3 -std=c99 -march=native

TARGETS  = pi_seq pi_omp
all: $(TARGETS) 

pi_seq: pi_seq.c
	$(CC) $(CFLAGS) $< -o $@

pi_omp: pi_omp.c
	$(CC) $(CFLAGS) $(OPENMP) $< -o $@

clean:
	rm -rf $(TARGETS)

ultraclean:
	rm -rf $(TARGETS) *.o* *.e* *.txt *.ps *.pdf time-*
