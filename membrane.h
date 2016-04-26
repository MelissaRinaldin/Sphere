#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef TYPES_DECLARATION

typedef struct {

  long d;
  long h;
  long m;
  long s;

} CPU_Time;

#endif
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define TYPES_DECLARATION
#define PI 3.141592653589793238462643383279502884
#define SQRT3 1.732050807568877293527446341505872366
