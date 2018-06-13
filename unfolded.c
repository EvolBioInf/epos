/***** unfolded.c *********************************
 * Description: Implement Peter Pfaffelhuber's
 *   optimization approach to population size esti-
 *   mation based on the unfolded SFS as described
 *   in his memo dated December 19, 2017.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Dec 18 09:51:48 2017
 **************************************************/
#include <stdio.h>
#include <float.h>
#include <gsl/gsl_linalg.h>
#include "util.h"
#include "sfs.h"
#include "popSizes.h"
#include "printMatrix.h"
#include "eprintf.h"

/* getCoeffMat: Left-hand side of equation (10) */
gsl_matrix *unfoldedGetCoeffMat(Sfs *sfs, PopSizes *ps, Args *args){
  gsl_matrix *mat;
  int i, j, n, m, *k, r;
  double a, factNki, u;
  double lambda;

  lambda = sfs->l;
  u = sfs->u;

  n = sfs->n;
  m = ps->m;
  k = ps->k;
  mat = gsl_matrix_alloc(m, m);
  for(j=1; j<=m; j++){
    for(i=1; i<=m; i++){
      factNki = 0.;
      for(r=1; r<= n-1; r++){
	a = binomial(n-1, r);
	factNki += 4.*u/r * 1./a/a * fi(j, n, r, k) * fi(i, n, r, k);
      }
      factNki += lambda * ((2-delta(j,1) - delta(j,m)) * delta(i,j) - delta(i-1,j) - delta(i+1,j));
      gsl_matrix_set(mat, j-1, i-1, factNki);
    }
  }
  return mat;
}

/* getResVect: Right-hand side of equation (10) */
gsl_vector *unfoldedGetResVect(Sfs *sfs, PopSizes *ps){
  gsl_vector *b;
  double x, a;
  int j, m, n, r, *k;

  b = gsl_vector_alloc(ps->m);
  n = sfs->n;
  m = ps->m;
  k = ps->k;
  for(j=1; j<=m; j++){
    x = 0.;
    for(r=1; r<=n-1; r++){
      a = 1. / binomial(n-1, r);
      x += sfs->f[r-1] * a * fi(j, n, r, k);
    }
    gsl_vector_set(b, j-1, x);
  }
  return b;
}

/* psi: First equation in Section 2.3 */
double unfoldedPsi(PopSizes *ps, Sfs *sfs){
  int i, n, m, r;
  double a, s1, s2, *N, u, lambda;

  s1 = 0.;
  m = ps->m;
  N = ps->N;
  n = sfs->n;
  u = sfs->u;
  lambda = sfs->l;
  for(r=1; r<=n-1; r++){
    s2 = 0.;
    for(i=1; i<=m; i++)
      s2 += N[i-1] * fi(i, n, r, ps->k);
    a = sfs->f[r-1] - 4.*u/r * 1./binomial(n-1,r) * s2;
    s1 += r * a * a;
  }
  s2 = 0.;
  a = 1;
  for(i=1; i<m; i++){
    a = N[i] - N[i-1];
    a *= a;
    s2 += a;
  }
  a = exp(log(a) + log(lambda));
  s1 += a;

  if(s1 < 0){
    printf("ERROR[unfoldedPsi]: Psi < 0: %g\n", s1);
    exit(-1);
  }

  return s1;
}
