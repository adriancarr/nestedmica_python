#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define BITS 12

#ifdef NMICA_BIG_ENDIAN
#define WORD0 0
#define WORD1 1
#define ENDIAN_COMMENT "big-endian"
#else
#define WORD0 1
#define WORD1 0
#define ENDIAN_COMMENT "little-endian"
#endif

static double* logtable;

void init_fastlog() {
	int top;
    volatile double x;
    volatile double *buffer = &x;
    volatile int *istar = (int *)buffer;
    istar[WORD1] = 0;
    double vlog2 = log(2.0);
    
    /* fprintf(stderr, "Initializing fastlog code.  Compiled for a %s platform.  Precision=%d.\n", ENDIAN_COMMENT, BITS); */
    
    logtable = malloc((1 << BITS) * sizeof(double));
    
    for (top = 0; top < (1 << BITS); ++top) {
        istar[WORD0] = (1023 << 20) | (top << 20 - BITS);
        logtable[top] = log(x) / vlog2;
    }
}

double fastlog2(volatile double f) {
    volatile int *istar = (int *) &f;
    int i0 = istar[WORD0];
    int e = (i0 >> 20) & ((1 << 11) - 1);
    int tf = (i0 >> (20 - BITS)) & ((1 << BITS) - 1);
    return logtable[tf] + e - 1023;
}
