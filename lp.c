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

#include "mt19937p.h"

/* References: [1] http://arxiv.org/abs/1407.1925 */

/* Currently getting bad results with float - seems like we might be dealing with values that are too small */
/* #define real float */
/* #define real_max FLT_MAX */
/* #define real_eps FLT_EPSILON */

#define real double
#define real_max DBL_MAX
#define real_eps DBL_EPSILON

#define AA(i,j) a[i +j*m]
/* b is a length-m vector */
/* c is a length-n vector*/
/* a is an mxn matrix in column major format */

real thresh(real v, real eps)
{
  /* This should be true, but this is wayyy in the inner loop */
  /* assert(v >= -1); */
  if (v > 1)
    {
      return 1;
    }
  else if ( v >= -eps && v <= eps)
    {
      return 0;
    }
  else
    return v;
}

void print_vec( real const * const x, int const n, int const stride)
{
  for(int idx = 0; idx < n-1; ++idx)
    printf("%f, ", x[idx*stride]);
  printf("%f\n", x[(n-1)*stride]);
}

void print_mat(real  const * const a, int const m, int const n)
{
  for (int ii = 0; ii < m; ++ii)
    {
      print_vec(a+ii, n, m);
    }
}

/*  Generates a random A matrix where entries are nonzero w/p p. */
/*  Except for first row, which has uniform [0,1] entries */
void random_instance(int const m, int const n, real ** a,real const p)
{

  /* int *l = (int*) _mm_malloc(n*n*sizeof(int), ALIGNBY); */
  real *l = (real*) malloc(m*n*sizeof(real));
  struct mt19937p state;
  sgenrand(10302011UL, &state);
  for (int i = 0; i < m; ++i){
    for (int j = 0; j < n; ++j) 
      {
	if( i == 0)
	  l[j*m+i] = genrand(&state);
	else 
	  l[j*m+i] = (genrand(&state) < p);
      }
  }
  /* return l; */
  *a = l;
}

int count_nonzeros(real const * const x, int const n)
{
  int count = 0;

  for(int ii = 0; ii < n; ++ii)
    {
      if( fabs(x[ii]) > real_eps)
	++count;
    }
  return count;
}

/* Initializes variables to specified size */
void gen_test_problem(real** A, real ** B, 
		      real ** C, int const m, int const n)
{
  real * a;
  real * b;
  real * c;

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

  *A = a;
  *B = b;
  *C = c;
}
/* Convert to problem of the form, say max 1^{T}x s.t. Ax \leq{} 1 */
void to_standard_form(real * const a, real const *const b, real const * const c, 
		      int const m, int const n)
{
  /* Scale down by b */
  for(int ii = 0; ii < m; ++ii)
    {
      for (int jj = 0; jj < n; ++jj)
	{
	  AA(ii,jj) = AA(ii,jj) / b[ii];
	}
    }
  
  /* Scale down by c */
  for(int ii = 0; ii < m; ++ii)
    {
      for (int jj = 0; jj < n; ++jj)
	{
	  AA(ii,jj) = AA(ii,jj) / c[jj];
	}
    }
}

/* Convert back to  max c^{T}x s.t. Ax \leq{} b */
void from_standard_form(real * const a, real const *const b, real const * const c,
			int const m, int const n)
{
  /* Scale up by b */
  for(int ii = 0; ii < m; ++ii)
    {
      for (int jj = 0; jj < n; ++jj)
	{
	  AA(ii,jj) = AA(ii,jj) * b[ii];
	}
    }
  
  /* Scale up by c */
  for(int ii = 0; ii < m; ++ii)
    {
      for (int jj = 0; jj < n; ++jj)
	{
	  AA(ii,jj) = AA(ii,jj) * c[jj];
	}
    }
}

real infinity_norm(real const * const x, int const n, int const stride)
{
  assert (stride > 0);

  real lb = 0;

  for(int idx = 0; idx < n; ++idx)
    {
      real const v = fabs(x[idx*stride]);
      if(v>lb)
	lb = v;
    }
  return lb;
}

int  infinity_norm_argmax(real const * const x, int const n, int const stride)
{
  assert (stride > 0);

  real lb = 0;

  int arg;

  for(int idx = 0; idx < n; ++idx)
    {
      real const v = fabs(x[idx*stride]);
      if(v>lb)
	{
	  lb = v;
	  arg = idx;
	}
    }
  return arg;
}

/* Compute Ax given A and x, where A is mxn and x is length-n*/
void matvec_multr(real const * const x, real const * const a, real * const ax, int const m, int const n)
{
  for(int ii = 0; ii < m; ++ii)
    {
      ax[ii] = 0.0;
    }

  for(int jj = 0; jj < n; ++jj)
    {
      for(int ii = 0; ii < m; ++ii)
	{
	  ax[ii] += AA(ii,jj)*x[jj];
	}
    }
}

/* Compute x^{T}A given A and x, where A is mxn and x is length-m */
void matvec_multl(real const * const x, real const * const a, real * const xa, int const m, int const n)
{
  for(int ii = 0; ii < n; ++ii)
    {
      xa[ii] = 0.0;
    }

  for(int jj = 0; jj < n; ++jj)
    {
      for(int ii = 0; ii < m; ++ii)
	{
	  xa[jj] += AA(ii,jj)*x[ii];
	}
    }
}

/* Exponentiate vector inplace (i.e., x <- e^(scale*(x+offset))*/
void exp_vec(real * const x, real const scale, real const offset, int const n)
{
  for (int ii = 0; ii < n; ++ii)
    {
      x[ii] = exp(scale*(x[ii]+offset));
    }
}

/* x1 <- x1 + x2 */
void add_vec(real * const x1, real * const x2, int const n)
{
  for (int ii = 0; ii < n; ++ii)
    {
      x1[ii] +=x2[ii];
    }
}
void mult_vec(real * const x1, real * const x2, int const n)
{
  for (int ii = 0; ii < n; ++ii)
    {
      x1[ii] *=x2[ii];
    }
}
void scale_vec(real * const x, real const scale, int const n)
{
  for (int ii = 0; ii < n; ++ii)
    {
      x[ii] *=scale;
    }
}
/* a is assumed to be column-major */
void scale_mat(real * const a, real const scale, int const m, int const n)
{
  for (int ii = 0; ii < n; ++ii)
    {
      scale_vec(a+ii*m,scale, m);
    }
}

/* Compute dual variable y given x */
void dual_ao15(real const * const x, real * const y, real const * const a, 
	       real mu, int const m, int const n)
{
  /* Set y = Ax */
  matvec_multr(x,a,y,m,n);

  exp_vec(y, 1/mu, -1, m);
}


/* See [1] page 19 */
void fixcoord_ao15(real * const y, real const * const a, real const eps, int const m, int const n)
{
  real * ya = (real*) malloc(n*sizeof(real));
  matvec_multl(y,a,ya,m,n);
  for(int ii = 0; ii < n; ++ii)
    {
      real const lambda = ya[ii] -1 + eps;
      if( lambda < -eps)
	{
	  int const jj = infinity_norm_argmax(a+ii*m, m, 1);
	  y[jj] = y[jj] - lambda/AA(jj,ii);
	}
    }

  scale_vec(y, 1.0/(1-2*eps), m);

  free(ya);
  ya = NULL;
}

/* Get duality gap for x and fixed y */
real test_duality_gap_fixed(real const * const a, real const * const x, real const * const y_in,
		      real const eps, int const m, int const n)
{
  real * y = (real*) malloc(m*sizeof(real));
  memcpy(y, y_in, m*sizeof(real));

  fixcoord_ao15(y, a, eps, m, n);

  real p = 0, d = 0;

 for(int ii =0; ii < n; ++ii)
    {
      p += x[ii];
    }
  for(int ii =0; ii < m; ++ii)
    {
      d += y[ii];
    }
  real const gap = d/p;
  
  free(y);
  y = NULL;

  return gap;
}
/* Solve standard form LP */
/* Output should be (1-O(eps))-approximately optimal and satisfy constraints precisely */

void lpsolve_ao15(real ** X, real ** Y, real const epsi, real const * const a, int const m, int const n)
{
  real *x = NULL;
  real *y = NULL; 

  /* Alg. in paper doesn't solve to (1-eps) but rather (1-O(eps)). This guarantees our solution is (1-eps)-optimal */
  real const eps = epsi/6;	
  /* y scratch variable */
  real *yk = NULL;
  /* feedback variable */
  real * v = NULL;

  /* Solution variables */
  x = (real*) malloc(n*sizeof(real)); /* Primal variable */
  y = (real*) malloc(m*sizeof(real)); /* Dual variable */
  yk = (real*) malloc(m*sizeof(real)); /* Dual var. scratch */
  v = (real*) malloc(n*sizeof(real));

  /* Constants */
  real const mu = eps/(4*log(n*m/eps));
  real const alpha = eps*mu/4.0;
  long int const T = ((int)6*log(2*n)/(alpha*eps)+1);
  
  printf("Mu: %f, Alpha: %f, Eps: %f\n", mu, alpha, eps);
  printf("Running for %d iterations.\n", T);

  /* Initial value of x */
  for(int ii = 0; ii < n; ++ii)
    {
      real const ainf = infinity_norm(a+ii*m, m, 1);
      x[ii] = (1-eps/2)/(n*ainf);
    }
  /* Intialize cumulative sum of yk */
  for(int ii = 0; ii < m; ++ii)
    {
      y[ii] = 0;
    }

  double start = omp_get_wtime();

  /* Main loop */
  for(long int t = 0; t < T; ++t)
    {
      if( !(t % 5000000))
      /* if( !(t % 500000)) */
	{
	  double now = omp_get_wtime();

	  real gap = test_duality_gap_fixed(a, x, y, eps, m, n);
	  printf("Reached iteraton %d (%f%%). Elapsed time %f\n", t, (float)t/T*100, now-start);
	  printf("Gap: %f, Expected: %f\n", gap, (1+epsi)/(1-epsi));
	  printf("Primal: ");
	  print_vec(x, n, 1);
	  printf("Dual: ");
	  print_vec(y, m, 1);
	  /* print_vec(yk, m, 1); */
	  fflush(stdout);
	}

      /* Compute y_{k} from x_{k} */
      dual_ao15(x,yk,a, mu, m,n);

      /* Set v = A^{T}y_{k}-1 */
      for(int jj = 0; jj < n; ++jj)
	{
	  v[jj]=-1;
	  for(int ii = 0; ii < m; ++ii)
	    {
	      v[jj] += AA(ii,jj)*yk[ii];
	    }
	}
      /* Update x^{k+1} */
      for(int jj = 0; jj < n; ++jj)
	{
	  x[jj] *= exp(-alpha*thresh(v[jj],eps));
	}

      /* Update cumulative sum */
      /* scale_vec(yk, 1./((double)T),m); */
      scale_vec(y, ((real)t)/(t+1), m);
      scale_vec(yk, 1./((real)(t+1)), m);
      add_vec(y,yk,m);
    }

  fixcoord_ao15(y, a, eps, m, n);
  /* At this point, y should be (1+4eps)OPT, so at most (1+epsi)OPT (since we scaled eps down) */

  scale_vec(x, 1/(1+eps),n);

  free(v);
  v = NULL;
  free(yk);
  yk = NULL;

  *X = x;
  *Y = y;
}

real min_colinf(real const *const a, int const m, int const n)
{
  real ub = real_max;
  for(int ii = 0; ii < n; ++ii)
    {
      real const ainf = infinity_norm(a+ii*m, m, 1);
      if (ub >= ainf)
	ub = ainf;
    }
  return ub;
}

void certify_standard_form(real const* const a, real const * const x, real const * const y, real const eps,
			   int const m, int const n)
{
  real p = 0, d=0;

  real * ax = (real*) malloc(m*sizeof(real));
  real * ya = (real*) malloc(n*sizeof(real));
  matvec_multr(x, a, ax, m,n);
  matvec_multl(y, a, ya, m,n);

  for(int ii = 0; ii < m; ++ii)
    {
      if (ax[ii] > 1 + real_eps)
	{
	  printf("Primal solution violates constraint %d (Value %f)\n", ii, ax[ii]);
	}
    }
  for(int ii = 0; ii <n ; ++ii)
    {
      if (ya[ii] < 1 - real_eps)
	{
	  printf("Dual solution violates constraint %d (Value %f)\n", ii, ya[ii]);
	}
    }
  for(int ii =0; ii < n; ++ii)
    {
      p += x[ii];
    }
  for(int ii =0; ii < m; ++ii)
    {
      d += y[ii];
    }

  printf("Primal value: %f, Dual value %f\n", p, d);

  real const gap = d/p;
  real const ub = (1+eps)/(1-eps);

  printf("Duality gap: %f, Expected: %f.\n",gap,ub);

  if( gap > ub)
    printf("Duality gap is too large.\n");

  free(ax);
  ax = NULL;
  free(ya);
  ya = NULL;
}

/* Solve packing LPs of the form max c^{T}x subject to Ax <= b */
int main(int argc, char **argv)
{
  int n = 100;			/* Dimensionality */
  int m = 200;			/* Number of constraints */
  /* int n = 2; */
  /* int m = 2; */
  int N = 0;			/* Number of non-zeros */
  real eps = 0.099;		/* precision */
  real *a = NULL, *b = NULL, *c = NULL;		/* instance */
  real *x = NULL, *y = NULL; 			/* solution variables */
  

  /* Read in parameters */


  printf("Variables: %d, Constraints: %d\n", n, m);

  random_instance(m,n, &a, 0.6);
  N = count_nonzeros(a,m*n);
  /* printf("%d non-zeros in A.\n", N); */

  /* Get test instance and convert to standard LP formulation */
  /* gen_test_problem(&a,&b,&c,m,n); */
  /* Put in standard packing LP form */
  /* to_standard_form(a,b,c,m,n); */
  /* This will totally break if b and c have zeros. Have to do a different reduction in this case */

  printf("A:\n");
  print_mat(a,m,n);

  /* Scale A down by infinity norm (as per page 5)  */
  /* Assumes A has no all-zero columns */
  real ainf = min_colinf(a, m,n);
  scale_mat(a, 1/ainf, m, n);

  /* Solve LP */
  assert(eps > 0 && eps <= 0.1);
  
  printf("Solving to precision %f.\n", eps);

  lpsolve_ao15(&x,&y,eps,a,m,n);
  /* Expect x and y to satisfy constraints perfectly and to have */
  /* 1^{T}x >= (1-eps)OPT, 1^{T}y <= (1+eps)OPT */
  printf("Primal solution:\n");
  print_vec(x,n,1);
  printf("Dual solution:\n");
  print_vec(y,m,1);

  /* /\* Maybe certify x using y here *\/ */
  certify_standard_form(a, x, y, eps, m, n);

  /* /\* Scale up solution (undo standardization and ainf scaling) *\/ */
  /* /\* for(int ii = 0; ii < n; ++ii) *\/ */
  /* /\*   { *\/ */
  /* /\*     x[ii] *= c[ii]*ainf; *\/ */
  /* /\*   } *\/ */

  free(x);
  x = NULL;
  free(y);
  y = NULL;
  free(a);
  a = NULL;
  /* free(b); */
  /* b = NULL; */
  /* free(c); */
  /* c = NULL; */
  return 0;
}
