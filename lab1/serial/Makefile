CC = gcc
OPENMP	 = 
CFLAGS   = -Wall -O3 -std=c99 -march=native

TARGETS  = pi_seq
all: $(TARGETS) 

pi_seq: pi_seq.c
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf $(TARGETS)

ultraclean:
	rm -rf TRACE*.mpit $(TARGETS) *.prv *.pcf *.row set-0 *.mpits *.sym *.o* *.e* time-* *.ps *.pdf *.fig *.txt
