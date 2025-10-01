all: BackSubs_mpi BackSubs_mpi_omp

BackSubs_mpi: BackSubs_mpi_omp.c
	mpicc -o BackSubs_mpi BackSubs_mpi_omp.c

BackSubs_mpi_omp: BackSubs_mpi_omp.c
	mpicc -fopenmp -o BackSubs_mpi_omp BackSubs_mpi_omp.c

clean:
	rm BackSubs_mpi BackSubs_mpi_omp