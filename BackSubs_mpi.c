/* C Example of Backward Substitution */
#define _XOPEN_SOURCE
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* header files for getting hostname and process id */
#include <sys/types.h>
#include <unistd.h>

/* Get time and resources */
#include <sys/resource.h>
#include <sys/time.h>

// #define _DEBUG_
// #include "constants.h"

#if _EXTRAE_
#include "extrae_user_events.h"
// Extrae Constants
#define PROGRAM 1000
#define END 0
#define SERIAL 1
#define PARALLEL 2
#else
double getusec_() {
    struct timeval time;
    gettimeofday(&time, NULL);
    return ((double)time.tv_sec * (double)1e6 + (double)time.tv_usec);
}

#define START_COUNT_TIME stamp = getusec_();
#define STOP_COUNT_TIME(_m)     \
    stamp = getusec_() - stamp; \
    stamp = stamp / 1e6;        \
    printf("Time: %s%0.6f\n", (_m), stamp);
#endif

#ifdef SINGLE_PRECISION

typedef float fp_t;

#else

typedef double fp_t;

#endif

int m;     /* order of symmetric input matrix */
int n;     /* number of right-hand sides */
int b;     /* block size */
int reps;  /* number of repetitions */
int check; /* check result */
fp_t alpha;
int lda;
int ldb;
char uplo; /* lower/upper triangular */
char side;
char trans;
char diag;

fp_t *A;
fp_t *B;
fp_t *Bchk;

#define SIZE 4

// Function to print the 4x4 matrix A and the 4-element vector B side-by-side
void print_matrix_vector(double A[SIZE][SIZE], double B[SIZE]) {
    printf("--- 4x4 Matrix (A) | 4-Element Vector (B) ---\n");

    // Outer loop for rows (0 to 3)
    for (int i = 0; i < SIZE; i++) {
        // 1. Print Matrix A row
        for (int j = 0; j < SIZE; j++) {
            // %8.3f for formatted output (total width 8, 3 decimal places)
            printf("%8.3f ", A[i][j]);
        }

        // 2. Print Separator and Vector B element
        printf("| %8.3f\n", B[i]);
    }
}

// Function to print the 4x4 matrix A and the 4-element vector B side-by-side
void print_matrix_vector_fake(int id, double A[2][SIZE], double B[2]) {
    printf("--- id: %d | 2x4 Matrix (A) | 2-Element Vector (B) ---\n", id);

    // Outer loop for rows (0 to 3)
    for (int i = 0; i < 2; i++) {
        // 1. Print Matrix A row
        for (int j = 0; j < SIZE; j++) {
            // %8.3f for formatted output (total width 8, 3 decimal places)
            printf("%8.3f ", A[i][j]);
        }

        // 2. Print Separator and Vector B element
        printf("| %8.3f\n", B[i]);
    }
}

void GENMAT_IP(fp_t *A, int m, int n, int scale, int seed) {
    srand48(time(NULL) + seed);

    int j;
    for (j = 0; j < n; ++j) {
        int i;
        for (i = 0; i < m; ++i) {
            A[j * m + i] = (fp_t)drand48();
            if (j == i) A[j * m + i] += scale;
        }
    }
}

int trsm_setup(int check, int m, int n, int b, int lda, int ldb, fp_t **A, fp_t **B, fp_t **Bchk) {
    fp_t *lA = *A = malloc(lda * m * sizeof(fp_t));
    if (lA == NULL) return 1;

    GENMAT_IP(lA, lda, m, 1, 0);
    int j, i;
    for (j = 0; j < m; ++j)
        for (i = 0; i < j; ++i) lA[j * m + i] = (fp_t)0.0;

    fp_t *lB = *B = malloc(ldb * n * sizeof(fp_t));
    if (lB == NULL) return 2;
    GENMAT_IP(lB, ldb, n, 1, 1);

    fp_t *lBchk = *Bchk = malloc(ldb * n * sizeof(fp_t));
    if (lBchk == NULL) return 3;

    for (i = 0; i < ldb * n; ++i) lBchk[i] = lB[i];

    return 0;
}

void trsm_shutdown(fp_t *A, fp_t *B, fp_t *X) {
    free(A);
    free(B);
    free(X);
}

int main(int argc, char *argv[]) {
    int myid, numprocs;

    MPI_Status status;

    const char Usage[] = "Usage: BackSubs <size> (try 10000)\n";
    if (argc < 2) {
        fprintf(stderr, Usage);
        exit(1);
    }
    unsigned long long int size = atoll(argv[1]);

    side = 'L';
    uplo = 'U';

    m = size;
    n = size;
    b = 1;
    lda = size;
    ldb = 1;

    alpha = 1.0;
    trans = 'N';
    diag = 'N';
    reps = 1;
    check = 1;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    int FACTOR = 2;

    if (numprocs < 2 || size % (FACTOR * numprocs)) MPI_Abort(MPI_COMM_WORLD, 1);

    if (myid == 0) {  // Only Id 0 creates the matrix to share it
        if (trsm_setup(check, m, n, b, lda, ldb, &A, &B, &Bchk)) {
            fprintf(stderr, "err: allocating matrix\n");
            return 2;
        }
    }

    // distribuzione carico
    // TODO non gestisco il resto

    int factor_processes = FACTOR * numprocs;
    int lines = size / factor_processes;
    int block_size = lines * size;

    fp_t *mem_A = malloc(FACTOR * block_size * sizeof(fp_t));
    fp_t *mem_B = malloc(FACTOR * lines * sizeof(fp_t));

    if (myid) {
        for (int i = FACTOR - 1; i >= 0; --i) {
            MPI_Recv(&mem_A[i * block_size], block_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&mem_B[i * lines], lines, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        }
    } else {
        for (int i = size - lines, recipient = 1; i >= 0; i -= lines, recipient++) {
            int dest = recipient % numprocs;
            if (dest != 0) {
                MPI_Send(&A[i * size], block_size, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
                MPI_Send(&B[i], lines, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
            } else {
                int j = FACTOR - recipient / numprocs;
                memcpy(&mem_A[j * block_size], &A[i * size], block_size * sizeof(fp_t));
                memcpy(&mem_B[j * lines], &B[i], lines * sizeof(fp_t));
            }
        }
    }

    // a questo punto tutti i processi hanno
    // mem_A = (FACTOR * block_size) coefficienti
    // mem_B = (FACTOR * lines) controvalori

    // time
    double stamp;
    if (myid == 0) START_COUNT_TIME;

    // la computazione effettiva comincia da proc 1 così finisce in proc 0

    int actual_proc = 1;
    unsigned global_index = size - lines;
    unsigned segment_to_be_updated = FACTOR;
    fp_t mem;

    for (int i = FACTOR - 1; i >= 0; --i) {
        for (int p = 0; p < numprocs; ++p) {
            if (myid == actual_proc) {
                unsigned offset_i = i * lines;
                for (int j = lines - 1; j >= 0; --j) {
                    mem_B[offset_i + j] = mem_B[offset_i + j] / mem_A[i * block_size + j * size + global_index + j];
                    MPI_Bcast(&mem_B[offset_i + j], 1, MPI_DOUBLE, actual_proc, MPI_COMM_WORLD);
                    for (int k = offset_i + j - 1; k >= 0; --k) {
                        mem_B[k] -= mem_B[offset_i + j] * mem_A[k * size + global_index + j];
                    }
                    if (myid == 0) B[global_index + j] = mem_B[offset_i + j];
                }
                segment_to_be_updated--;
            } else {
                for (int j = lines - 1; j >= 0; --j) {
                    MPI_Bcast(&mem, 1, MPI_DOUBLE, actual_proc, MPI_COMM_WORLD);
                    for (int k = segment_to_be_updated * lines - 1; k >= 0; --k) {
                        mem_B[k] -= mem * mem_A[k * size + global_index + j];
                    }
                    if (myid == 0) B[global_index + j] = mem;
                }
            }
            actual_proc = (actual_proc + 1) % numprocs;
            global_index -= lines;
        }
    }

    if (myid == 0) {
        STOP_COUNT_TIME("");
    }

    if (myid == 0) {
        if (size <= 1000) {
            if (size <= 10) {  // Simple printing code for debug purposes
                for (int i = 0; i < lda; ++i) {
                    for (int j = 0; j < m; ++j) fprintf(stdout, "%lf ", A[i * m + j]);
                    fprintf(stdout, "\n");
                }
                fprintf(stdout, "\n");
                for (int i = 0; i < ldb * n; ++i) fprintf(stdout, "%lf ", B[i]);
                fprintf(stdout, "\n");
                for (int i = 0; i < ldb * n; ++i) fprintf(stdout, "%lf ", Bchk[i]);
                fprintf(stdout, "\n");
            }
            fp_t error = 0;
            fp_t suma;
            for (int i = size - 1; i >= 0; i--) {
                suma = 0.0;
                for (int j = i; j <= size - 1; j++) suma += B[j] * A[i * m + j];
                if (error < fabs(Bchk[i] - suma)) {
                    error = fabs(Bchk[i] - suma);
                    fprintf(stdout, "%d:%lf vs %lf -> %lf\n", i, suma, Bchk[i], error);
                }
            }
            fprintf(stdout, "Error: %lf\n", error);
        }
    }

    if (myid == 0) {
        trsm_shutdown(A, B, Bchk);
    }
    // free
    free(mem_A);
    free(mem_B);

    MPI_Finalize();

    return EXIT_SUCCESS;
}
