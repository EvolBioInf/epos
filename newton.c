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

Rparams *newRparams(Sfs *sfs, PopSizes *ps) {
  Rparams *r = (Rparams *)emalloc(sizeof(Rparams));

  r->n = sfs->n;
  r->m = ps->m;
  r->k = ps->k;
  r->g = sfs->f;
  r->u = sfs->u;

  return r;
}

/* expG implements equation (2) */
double expG(double *N, double u, int m, int n, int *k, int r) {
  int i;
  double a, b, s, x;

  s = 0.;
  for(i = 0; i < m; i++) {
    x = n - k[i] + 1;
    if(x < 0)
      x = 0;
    a = binomial(x, r);
    x = n - k[i+1] + 1;
    if(x < 0)
      x = 0;
    b = binomial(x, r);
    s += N[i] * (a - b);
  }
  s *= 4. * u / (double)r / binomial(n - 1, r);
  
  return s;
}

int func(const gsl_vector *x, void *params, gsl_vector *f) {
  int i, r;
  double eg, xx, yy;
  int n     = ((Rparams *) params)->n;
  int m     = ((Rparams *) params)->m;
  int *k    = ((Rparams *) params)->k;
  double *g = ((Rparams *) params)->g;
  double u  = ((Rparams *) params)->u;

  double *N = (double *)emalloc(m * sizeof(double));
  double *y = (double *)emalloc(n * sizeof(double));
  
  for(i = 0; i < m; i++) {
    N[i] = gsl_vector_get(x, i);
  }

  for(i = 0; i < m; i++) {
    y[i] = 0.;
    for(r = 1; r < n; r++) {
      /* equation (3) */
      eg = expG(N, u, m, n, k, r);
      xx = n - k[i] + 1;
      if(xx < 0)
	xx = 0;
      yy = n - k[i+1] - 1;
      if(yy < 0)
	yy = 0;
      y[i] += 1. / r * (g[r-1] / eg - 1) * (binomial(xx, r) - binomial(yy, r)) / binomial(n - 1, r);
    }
  }
  for(i = 0; i < m; i++) {
    gsl_vector_set(f, i, y[i]);
  }

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

double logLik(PopSizes *ps, Sfs *sfs) {
  double e, l;

  l = 0.;
  for(int r = 1; r < ps->n; r++) {
    e = expG(ps->N, sfs->u, ps->m, ps->n, ps->k, r);
    l += sfs->f[r-1] * log(e) - e;
  }

  return l;
}

int newtonComp(Sfs *sfs, PopSizes *ps, Args *args, double iniP) {
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;
  int status;
  size_t i, iter = 0;
  Rparams *p = newRparams(sfs, ps);
  gsl_multiroot_function f = {&func, ps->m, p};
  gsl_vector *x = gsl_vector_alloc(ps->m);

  for(i = 0; i < ps->m; i++) {
    gsl_vector_set(x, i, iniP);
  }
  T = gsl_multiroot_fsolver_hybrids;
  s = gsl_multiroot_fsolver_alloc(T, ps->m);
  gsl_multiroot_fsolver_set(s, &f, x);

  status = 0;
  do {
    iter++;
    status = gsl_multiroot_fsolver_iterate(s);
    if(status) {
      fprintf(stderr, "# WARNING: The solver is stuck at iteration %ld.\n# Status of solver: ", iter);
      printState(iter, s, ps->m);
      return status;
    }
    status = gsl_multiroot_test_residual(s->f, 1e-7);
  } while (status == GSL_CONTINUE && iter < 1000);
  for(i = 0; i < ps->m; i++) {
    ps->N[i] = gsl_vector_get(s->x, i);
  }
  gsl_multiroot_fsolver_free(s);
  gsl_vector_free(x);

  return status;
}

int newton(Sfs *sfs, PopSizes *ps, Args *args) {
  int status;
  double iniP, origP;

  iniP = watterson(sfs);
  origP = iniP;
  status = newtonComp(sfs, ps, args, iniP);

  while(status && iniP > origP / 10e4) {
    iniP /= 2.;
    fprintf(stderr, "# Trying again with %f as initial population size.\n", iniP);
    status = newtonComp(sfs, ps, args, iniP);
  }

  if(status) {
    fprintf(stderr, "# ERROR: Solver failed, aborting computation.\n");
    exit(-1);
  }

  return status;
}
