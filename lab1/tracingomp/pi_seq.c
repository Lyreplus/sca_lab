/*
 * Compute pi by Monte Carlo
 */
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#if _EXTRAE_
#include "extrae_user_events.h"
// Extrae Constants
#define  PROGRAM    1000
#define  END        0
#define  SERIAL     1
#define  PARALLEL   2
#else 
double getusec_() {
        struct timeval time;
        gettimeofday(&time, NULL);
        return ((double)time.tv_sec * (double)1e6 + (double)time.tv_usec);
}


#define START_COUNT_TIME stamp = getusec_();
#define STOP_COUNT_TIME(_m) stamp = getusec_() - stamp;\
                        stamp = stamp/1e6;\
                        printf ("%s%0.6f\n",(_m), stamp);
#endif

inline float my_rand(unsigned long long int *seed) { 
       unsigned long long int a = 16807;  // constants for random number generator 
       unsigned long long int m = 2147483647;   // 2^31 - 1 
       unsigned long long int x = (unsigned long long int ) *seed; 
       x = (a * x)%m; 
       *seed = (unsigned long long int) x; 
       return ((float)x)/m; 
} 

int main(int argc, char *argv[]) {
#if _EXTRAE_
    Extrae_event (PROGRAM, SERIAL);
#else
    double stamp;
    START_COUNT_TIME;
#endif

    float x, y, pi=0.0;
    unsigned long long int points_in_circle=0;

    const char Usage[] = "Usage: pi <trials> (try 1000000000)\n";
    if (argc < 2) {
	fprintf(stderr, Usage);
	exit(1);
    }
    unsigned long long int trials = atoll(argv[1]);

#if _EXTRAE_
    Extrae_event (PROGRAM, END);
#endif

    /* do computation -- using all available threads */
#if _EXTRAE_
    Extrae_event (PROGRAM, PARALLEL);
#endif
     unsigned long long int seed = 0 + 1;
     for(unsigned long long int i = 0; i < trials; i++) {
	x = my_rand(&seed);
	y = my_rand(&seed);
	points_in_circle += (x*x + y*y <= 1.0f);
     }
#if _EXTRAE_
    Extrae_event (PROGRAM, END);
    Extrae_event (PROGRAM, SERIAL);
#endif
    pi = 4.0f * points_in_circle / trials;

    /* print results */
    printf("Number pi after %lld iterations = %.15f\n", trials, pi);

#if _EXTRAE_
    Extrae_event (PROGRAM, END);
#else
    STOP_COUNT_TIME("");
#endif

    return EXIT_SUCCESS;
}
