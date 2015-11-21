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

/* References: [1] http://arxiv.org/abs/1407.1925 */

#define real float
#define real_max FLT_MAX
#define real_eps FLT_EPSILON

#define AA(i,j) a[i +j*m]
/* b is a length-m vector */
/* c is a length-n vector*
/* a is an mxn matrix in column major format */

real thresh(real v, real eps)
{
  assert(v >= -1);
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
  for(int ii = 0; ii < m; ++ii)
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
      x[ii] *scale;
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



/* Solve standard form LP */
/* Output should be (1-O(eps))-approximately optimal and satisfy constraints precisely */

void lpsolve_ao15(real ** X, real ** Y, real const epsi, real const * const a, int const m, int const n)
{
  real *x;
  real *y; 

  /* Alg. in paper doesn't solve to (1-eps) but rather (1-O(eps)). This guarantees our solution is (1-eps)-optimal */
  real const eps = epsi/6;	
  /* y scratch variable */
  real *yk;
  /* feedback variable */
  real * v;

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
  
  /* Main loop */
  for(long int t = 0; t < T; ++t)
    {
      if( !(t % 5000000))
	{
	  printf("Reached iteraton %d (%f%%)\n", t, (float)t/T*100);
	  print_vec(x, n, 1);
	}
      
      /* Compute y_{k} from x_{k} */
      dual_ao15(x,yk,a, mu, m,n);
      /* Update cumulative sum */
      add_vec(y,yk,m);

      /* Set v = A^{T}y_{k}-1 */
      for(int jj = 0; jj < n; ++jj)
	{
	  v[jj]=-1;
	  for(int ii = 0; ii < m; ++ii)
	    {
	      v[jj] += AA(ii,jj)*y[ii];
	    }
	}
      /* Update x^{k+1} */
      for(int jj = 0; jj < n; ++jj)
	{
	  x[jj] *= exp(-alpha*thresh(v[jj],eps));
	}
    }

  scale_vec(x, 1/(1+eps),n);

  free(v);
  free(yk);

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

  printf("Duality gap: %f, Expected: %f.\n",gap,ub);;

  if( gap > ub)
    printf("Duality gap is too large.\n");

  free(ax);
}

/* Solve packing LPs of the form max c^{T}x subject to Ax <= b */
int main(int argc, char **argv)
{
  int n = 5;			/* Dimensionality */
  int m = 5;			/* Number of constraints */
  real eps = 0.099;		/* precision */
  real *a, *b, *c;		/* instance */
  real *x, *y; 			/* solution variables */

  /* Read in parameters */

  /* Get test instance and convert to standard LP formulation */
  gen_test_problem(&a,&b,&c,m,n);


  printf("Variables: %d, Constraints: %d\n", n, m);
  /* Put in standard packing LP form */
  to_standard_form(a,b,c,m,n);

  printf("A:\n");
  print_mat(a,m,n);

  /* Scale A down by infinity norm (as per page 5)  */
  real ainf = min_colinf(a, m,n);
  scale_mat(a, 1/ainf, m, n);

  /* Solve LP */
  assert(eps > 0 && eps <= 0.1);
  
  printf("Solving to precision %f.\n", eps);

  lpsolve_ao15(&x,&y,eps,a,m,n);

  print_vec(x,n,1);
  print_vec(y,m,1);

  /* Maybe certify x using y here */
  certify_standard_form(a, x, y, eps, m, n);

  /* Scale up solution (undo standardization and ainf scaling) */
  for(int ii = 0; ii < n; ++ii)
    {
      x[ii] *= c[ii]*ainf;
    }

  free(x);
  free(y);
  free(a);
  free(b);
  free(c);
  return 0;
}
