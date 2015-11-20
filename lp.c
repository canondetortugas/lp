#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <omp.h>

#define real float

#define AA(i,j) a[i +j*m]
/* b is a length-m vector */
/* c is a length-n vector*
/* a is an mxn matrix in column major format */


/* Initializes variables to specified size */
void gen_test_problem(real* a, real * b, 
		      real * c, int const m, int const n)
{
  a = (real*) malloc(n*m*sizeof(real));
  b = (real*) malloc(m*sizeof(real));
  c = (real*) malloc(n*sizeof(real));

  for(int idx = 0; idx < n; ++idx)
    {
      for(int idy = 0; idy < m; ++idy)
	{
	  AA(idy,idx)=1;
	}
    }

  for(int idx = 0; idx < m; ++idx)
    {
      b[idx]=1;
    }

  for(int idx = 0; idx < n; ++idx)
    {
      c[idx]=1;
    }
}

int main(int argc, char **argv)
{
  int n = 10;
  int m = 100;

  real *a, *b, *c;

  gen_test_problem(a,b,c,m,n);
  

  free(a);
  free(b);
  free(c);
  return 0;
}
