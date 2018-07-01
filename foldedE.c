/***** foldedE.c **********************************
 * Description: Equations taken from Peter Pfaffel-
 *   huber's write-up dated December 19, 2017.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Dec 18 11:30:33 2017
 **************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <gsl/gsl_linalg.h>
#include "popSizes.h"
#include "util.h"
#include "sfs.h"
#include "printMatrix.h"
#include "eprintf.h"

/* getCoeffMat: Left-hand side of equation (12) */
gsl_matrix *foldedEgetCoeffMat(Sfs *sfs, PopSizes *ps, Args *args){
  gsl_matrix *mat;
  int i, j, n, m, *k, r;
  double a, b, u, lambda, factNki;

  n = sfs->n;
  m = ps->m;
  k = ps->k;
  u = sfs->u;
  lambda = sfs->l;
  mat = gsl_matrix_alloc(m, m);
  for(j=1; j<=m; j++){
    for(i=1; i<=m; i++){
      factNki = 0.;
      for(r=1; r<=n/2; r++){
	a = binomial(n-delta(2*r,n), r);
	b = binomial(n-1,r);
	factNki += 4.*u/(double)(r) * 1./a/b * hi(j, n, r, k) * hi(i, n, r, k);
      }
      factNki += lambda * ((2-delta(j,1) - delta(j,m)) * delta(i,j) - delta(i-1,j) - delta(i+1,j));
      gsl_matrix_set(mat, j-1, i-1, factNki);
    }
  }
  return mat;
}

/* getResVect: Right-hand side of equation (12) */
gsl_vector *foldedEgetResVect(Sfs *sfs, PopSizes *ps){
  gsl_vector *b;
  double x, a;
  int j, m, n, r, *k;

  b = gsl_vector_alloc(ps->m);
  n = sfs->n;
  m = ps->m;
  k = ps->k;
  for(j=1; j<=m; j++){
    x = 0.;
    for(r=1; r<=n/2; r++){
      a = 1. / binomial(n-delta(2*r,n), r);
      x += sfs->f[r-1] * a * hi(j, n, r, k);
    }
    gsl_vector_set(b, j-1, x);
  }
  return b;
}

/* foldedEpsi: Equation just above equation (10) on p. 7.
 *   The equation consists of the sum of three complex terms,
 *   which I call - from left to right - x, y, and z.
 */
double foldedEpsi(PopSizes *ps, Sfs *sfs){
  int i, n, m, r;
  double a, s1, s2, *N;
  double u, lambda, x, y, z, psi;

  m = ps->m;
  N = ps->N;
  n = sfs->n;
  u = sfs->u;
  lambda = sfs->l;

  /* calculate x */
  s1 = 0.;
  for(r=1; r<=n/2-1; r++){
    s2 = 0.;
    for(i=1; i<=m; i++)
      s2 += N[i-1] * gi(i, n, r, ps->k);
    a = sfs->f[r-1] - 4.*u/r * 1./binomial(n-1,r) * s2;
    s1 += (double)(r*(n-r)/n) * a * a;
  }
  x = s1;

  /* calculate y */
  s1 = 0.;
  for(i=1; i<=m; i++)
    s1 += N[i-1] * fi(i, n, n/2, ps->k);
  a = sfs->f[n/2-1] - 4.*u/(n/2.)*1./binomial(n-1,n/2) * s1;
  y = n/2. * a * a;

  /* calculate z */
  s1 = 0.;
  for(i=1; i<=m-1; i++){
    a = N[i] - N[i-1];
    s1 += a * a;
  }
  z = exp(log(s1) + log(lambda));

  psi = x + y + z;

  if(psi < 0){
    fprintf(stderr, "ERROR[foldedEpsi]: Psi < 0: %g\n", psi);
    exit(-1);
  }

  return psi;
}

double foldedEchiSquared(PopSizes *ps, Sfs *sfs){
  int i, n, m, r;
  double a, w, sw, xs, s1, s2, *N, varF;
  double u;

  m = ps->m;
  N = ps->N;
  n = sfs->n;
  u = sfs->u;

  s1 = 0.;
  sw = 0.;
  for(r=1; r<=n/2-1; r++){
    s2 = 0.;
    for(i=1; i<=m; i++)
      s2 += N[i-1] * gi(i, n, r, ps->k);
    a = sfs->f[r-1] - 4.*u/r * 1./binomial(n-1,r) * s2;
    varF = sfs->f[r-1] * (1 - sfs->f[r-1]) / n;
    w = (double) (r*(n-r)/n);
    s1 += w * a * a / varF;
    sw += w;
  }
  xs = s1 / sw;

  return xs;
}
