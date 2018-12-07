/***** newton.c ***********************************
 * Description: Use muti-dimensional root-finding
 *   (Newton's method) to solve equation (3) in
 *   Peter Pfaffelhuber's memo entitled "Some small
 *   remarks", which he sent on July 5th 2018.
 * Source: M. Galassi et al. (2005). GNU Scientific 
 *   Library Reference Manual, Edition 1.6, for GSL
 *   Version 1.6. p. 454f.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Tue Jul 24 11:42:40 2018
 **************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include "eprintf.h"
#include "util.h"
#include "sfs.h"
#include "popSizes.h"
#include "newton.h"
#include "tab.h"

Rparams *newRparams(Sfs *sfs, PopSizes *ps) {
  Rparams *r = (Rparams *)emalloc(sizeof(Rparams));
  r->s = sfs;
  r->p = ps;
  return r;
}

/* folded estimates the population sizes from a folded 
 * site frequency spectrum using equation (4b).
 */
int folded(const gsl_vector *x, void *params, gsl_vector *f) {
  double e, eg, bb;
  int a, b;

  Sfs      *s = ((Rparams *) params)->s;
  PopSizes *p = ((Rparams *) params)->p;
  int       n = s->n;
  int       m = p->m;
  int      *k = p->k;
  int      *G = s->G;
  double   *N = p->N;
  
  double *y = (double *)emalloc((m + 1) * sizeof(double));
  for(int i = 1; i <= m; i++)
    N[i] = gsl_vector_get(x, i - 1);
  e = expF(p, s, 0);
  for(int i = 1; i <= m; i++) {
    y[i] = 0.;
    a = max(n - k[i]   + 1, 0);
    b = max(n - k[i+1] - 1, 0);
    for(int r = 1; r <= n/2 - 1; r++) {
      eg = expF(p, s, r);
      bb = binomial(a, r) - binomial(b, r) + binomial(a, n-r) - binomial(b, n-r);
      y[i] += (G[r] / eg - G[0]/e) / (double)r * bb / binomial(n-1, r);
    }
    eg = expF(p, s, n/2);
    bb = binomial(a, n/2) - binomial(b, n/2);
    y[i] += 2./n * (G[n/2] / eg - G[0]/e) * bb / binomial(n-1, n/2);                                                                            
  }
  for(int i = 1; i <= m; i++)
    gsl_vector_set(f, i-1, y[i]);

  free(y);
  return GSL_SUCCESS;
}

/* unfolded computes population sizes for unfolded frequency spectra 
 * and returns the status indicator of the GSL method employed.
 * The mathematics is listed in equation (S5).
 */
int unfolded(const gsl_vector *x, void *params, gsl_vector *f) {
  double e, eg, bb;
  int a, b;

  Sfs      *s = ((Rparams *) params)->s;
  PopSizes *p = ((Rparams *) params)->p;
  int       n = s->n;
  int       m = p->m;
  int      *k = p->k;
  int      *G = s->G;
  double   *N = p->N;
  double   *y = (double *)emalloc((m + 1) * sizeof(double));
  
  for(int i = 1; i <= m; i++)
    N[i] = gsl_vector_get(x, i - 1);
  e = expG(p, s, 0);
  for(int i = 1; i <= m; i++) {
    y[i] = 0.;
    a = max(n - k[i]   + 1, 0);
    b = max(n - k[i+1] - 1, 0);
    for(int r = 1; r <= n-1; r++) {
      eg = expG(p, s, r);
      bb = binomial(a, r) - binomial(b, r);
      y[i] += 1. / r * (G[r] / eg - G[0] / e) * bb / binomial(n - 1, r);
    }
  }
  for(int i = 1; i <= m; i++)
    gsl_vector_set(f, i-1, y[i]);
  free(y);

  return GSL_SUCCESS;
}

void printState(int iter, gsl_multiroot_fsolver *s, int n) {
  int i;

  fprintf(stderr, "# iter: %d; x:", iter);
  for(i = 0; i < n; i++) {
    fprintf(stderr, " %e", gsl_vector_get(s->x, i));
  }
  fprintf(stderr, "; y:");
  for(i = 0; i < n; i++) {
    fprintf(stderr, " %e", gsl_vector_get(s->f, i));
  }
  fprintf(stderr, "\n");
}

int newton(Sfs *sfs, PopSizes *ps, Args *args) {
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;
  int status;
  size_t iter = 0;
  Rparams *p = newRparams(sfs, ps);
  gsl_multiroot_function f; 
  gsl_vector *x = gsl_vector_alloc(ps->m);

  if(sfs->f)
    f.f = &folded;
  else
    f.f = &unfolded;
  f.n = ps->m;
  f.params = p;
  double w = watterson(sfs);
  for(int i = 1; i <= ps->m; i++)
    gsl_vector_set(x, i-1, w);
  T = gsl_multiroot_fsolver_hybrids;
  s = gsl_multiroot_fsolver_alloc(T, ps->m);
  gsl_multiroot_fsolver_set(s, &f, x);
  status = 0;
  do {
    iter++;
    status = gsl_multiroot_fsolver_iterate(s);
    if(status) {
      gsl_multiroot_fsolver_free(s);
      gsl_vector_free(x);
      free(p);
      ps->l = logLik(ps, sfs);
      return status;
    }
    status = gsl_multiroot_test_residual(s->f, 1e-7);
  } while (status == GSL_CONTINUE && iter < 1000);
  for(int i = 1; i <= ps->m; i++) {
    ps->N[i] = gsl_vector_get(s->x, i - 1);
    if(ps->N[i] < 1)
      ps->N[i] = 1.; /* smallest population size possible: 1 individual */
  }
  ps->l = logLik(ps, sfs);
  gsl_multiroot_fsolver_free(s);
  gsl_vector_free(x);
  free(p);

  return status;
}

void testNewton() {
  char *fn = "data/kap144i.dat";
  FILE *fp = efopen(fn, "r");
  Args *args = newArgs();
  Sfs *sfs = readSfs(fp, args);
  fclose(fp);
  iniBinom(sfs->n);
  PopSizes *ps = newPopSizes(sfs);
  ps->m = 4;
  ps->k[1] = 2;
  ps->k[2] = 4;
  ps->k[3] = 6;
  ps->k[4] = 25;
  ps->k[ps->m + 1] = sfs->n + 1;
  newton(sfs, ps, args);
  freeArgs(args);
  printSfsStats(sfs);
  printTimes(ps, sfs);
  freeSfs(sfs);
  freePopSizes(ps);
  freeBinom();
  tabFree();
}
