#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <omp.h>

/* The crown jewel */
#include <mkl.h>

double m[] = {
  3, 1, 3,
  1, 5, 9,
  2, 6, 5
};

double x[] = {
  -1, -1, 1
};

double y[] = {
  0, 0, 0
};

int
main()
{
  int i, j;

  for (i=0; i<3; ++i) {
    for (j=0; j<3; ++j) printf("%5.1f", m[i*3+j]);
    putchar('\n');
  }

  cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, m, 3,
	      x, 1, 0.0, y, 1);

  for (i=0; i<3; ++i)  printf("%5.1f\n", y[i]);

  return 0;
}
