/*
 * Compute pi by Monte Carlo
 */
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

double getusec_() {
        struct timeval time;
        gettimeofday(&time, NULL);
        return ((double)time.tv_sec * (double)1e6 + (double)time.tv_usec);
}

inline float my_rand(unsigned long long int *seed) { 
       unsigned long long int a = 16807;  // constants for random number generator 
       unsigned long long int m = 2147483647;   // 2^31 - 1 
       unsigned long long int x = (unsigned long long int ) *seed; 
       x = (a * x)%m; 
       *seed = (unsigned long long int) x; 
       return ((float)x)/m; 
} 

#define START_COUNT_TIME stamp = getusec_();
#define STOP_COUNT_TIME(_m) stamp = getusec_() - stamp;\
                        stamp = stamp/1e6;\
                        printf ("%s%0.6f\n",(_m), stamp);

int main(int argc, char *argv[]) {
    double stamp;
    START_COUNT_TIME;

    float x, y, pi=0.0;
    unsigned long long int points_in_circle=0;

    const char Usage[] = "Usage: pi <trials> (try 1000000000)\n";
    if (argc < 2) {
	fprintf(stderr, Usage);
	exit(1);
    }
    unsigned long long int trials = atoi(argv[1]);


     unsigned long long int seed = 0+1;
     for(unsigned long long int i = 0; i < trials; i++) {
	x = my_rand(&seed);
	y = my_rand(&seed);
	points_in_circle += (x*x + y*y <= 1.0f);
     }
    pi = 4.0f * points_in_circle / trials;

    /* print results */
    printf("Number pi after %lld iterations = %.15f\n", trials, pi);

    STOP_COUNT_TIME("");

    return EXIT_SUCCESS;
}
