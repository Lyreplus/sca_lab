binaries=BackSubs_mpi BackSubs_mpi_omp

all: $(binaries)

BackSubs_mpi: BackSubs_mpi_omp.c
	mpicc -o BackSubs_mpi BackSubs_mpi_omp.c

BackSubs_mpi_omp: BackSubs_mpi_omp.c
	mpicc -fopenmp -o BackSubs_mpi_omp BackSubs_mpi_omp.c

BackSubs_mpi_omp_sanitizer: BackSubs_mpi_omp.c
	mpicc -fsanitize=address -fopenmp -o BackSubs_mpi_omp_sanitizer BackSubs_mpi_omp.c

clean:$(binaries)
	rm $(binaries)