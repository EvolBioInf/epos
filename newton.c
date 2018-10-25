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
  r->x = sfs->x;
  r->u = sfs->u;
  r->l = sfs->nullCount + sfs->numPol;
  r->o = sfs->nullCount;
  return r;
}

/* expGr implements equation (2) */
double expGr(double *N, double u, int m, int n, int l, int *k, int r) {
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
  s *= 4. * u * l / (double)r / binomial(n - 1, r);
  /* s *= 4. * u / (double)r / binomial(n - 1, r); */
  
  return s;
}

/* expG implements equation (1) */
double expG(double *N, double u, int m, int n, int l, int *k, int r){
  double x, s;
  
  if(r == 0){
    s = 0.;
    for(int i=1; i<n; i++){
      x = expGr(N, u, m, n, l, k, i);
      s += x;
    }
    x = (double)l - s;
    return x;
  }else{
    return expGr(N, u, m, n, l, k, r);
  }
}

double expF(double *N, double u, int m, int n, int l, int *k, int r) {
  double f;

  if(r == 0)
    f = expG(N, u, m, n, l, k, r);
  else
    f = expG(N, u, m, n, l, k, r) + expG(N, u, m, n, l, k, n-r);
  
  return f;
}

int folded(const gsl_vector *x, void *params, gsl_vector *f) {
  int i, r;
  double e, eg, xx, yy, bb, q;
  char excluded;
  int n     = ((Rparams *) params)->n;
  int m     = ((Rparams *) params)->m;
  int *k    = ((Rparams *) params)->k;
  double *g = ((Rparams *) params)->g;
  char *ex  = ((Rparams *) params)->x;
  double u  = ((Rparams *) params)->u;
  int l     = ((Rparams *) params)->l;
  int o     = ((Rparams *) params)->o;
  double *N = (double *)emalloc(m * sizeof(double));
  double *y = (double *)emalloc(n * sizeof(double));
  
  for(i = 0; i < m; i++) {
    N[i] = gsl_vector_get(x, i);
  }
  e = expF(N, u, m, n, l, k, 0);
  /* is a frequency category excluded? */
  excluded = 0;
  for(i = 1; i < n/2; i++)
    if(ex[i-1])
      excluded = 1;
  if(excluded)
    q = 1;
  else
    q = (double)o / e;
  /* equation (4b) */
  for(i = 0; i < m; i++) {
    y[i] = 0.;
    xx = n - k[i] + 1;
    if(xx < 0)
      xx = 0;
    yy = n - k[i+1] - 1;
    if(yy < 0)
      yy = 0;
    for(r = 1; r < n/2; r++) {
      if(ex[r-1]) /* is frequency category f[r-1] excluded from the analysis? */
	continue;
      eg = expF(N, u, m, n, l, k, r);
      bb = binomial(xx, r) - binomial(yy, r) + binomial(xx, n-r) - binomial(yy, n-r);
      y[i] += (g[r-1] / eg - q) / (double)r * bb / binomial(n - 1, r);
    }
    if(ex[n/2-1]) /* is frequency category f[n/2-1] excluded from the analysis? */
      continue; 
    eg = expF(N, u, m, n, l, k, n/2);
    bb = binomial(xx, n/2) - 2. * binomial(yy, n/2);
    y[i] += 2./n * (g[n/2 - 1] / eg - q) * bb / binomial(n-1, n/2);
  }
  for(i = 0; i < m; i++) {
    gsl_vector_set(f, i, y[i]);
  }

  free(N);
  free(y);
  return GSL_SUCCESS;
}

int unfolded(const gsl_vector *x, void *params, gsl_vector *f) {
  int i, r;
  char excluded;
  double e, eg, xx, yy, bb, q;
  int n     = ((Rparams *) params)->n;
  int m     = ((Rparams *) params)->m;
  int *k    = ((Rparams *) params)->k;
  double *g = ((Rparams *) params)->g;
  char *ex  = ((Rparams *) params)->x;
  double u  = ((Rparams *) params)->u;
  int l     = ((Rparams *) params)->l;
  int o     = ((Rparams *) params)->o;
  double *N = (double *)emalloc(m * sizeof(double));
  double *y = (double *)emalloc(n * sizeof(double));
  
  for(i = 0; i < m; i++) {
    N[i] = gsl_vector_get(x, i);
  }
  e = expG(N, u, m, n, l, k, 0);
  /* is a frequency category excluded? */
  excluded = 0;
  for(i = 1; i < n; i++)
    if(ex[i-1])
      excluded = 1;
  if(excluded)
    q = 1;
  else
    q = (double)o / e;
  /* equation (3) */
  for(i = 0; i < m; i++) {
    y[i] = 0.;
    xx = n - k[i] + 1;
    if(xx < 0)
      xx = 0;
    yy = n - k[i+1] - 1;
    if(yy < 0)
      yy = 0;
    for(r = 1; r < n; r++) {
      if(ex[r-1]) /* is frequency category f[r-1] excluded from the analysis? */
	continue;
      eg = expG(N, u, m, n, l, k, r);
      bb = binomial(xx, r) - binomial(yy, r);
      y[i] += 1. / r * (g[r-1] / eg - q) * bb / binomial(n - 1, r);
    }
  }
  for(i = 0; i < m; i++) {
    gsl_vector_set(f, i, y[i]);
  }
  free(N);
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

double logLik(PopSizes *ps, Sfs *sfs) {
  double e, l, x;
  int max;

  if(sfs->type == UNFOLDED)
    max = ps->n - 1;
  else if(sfs->type == FOLDED_EVEN)
    max = ps->n / 2;
  else {
    fprintf(stderr, "ERROR[epos]: Can only deal with unfolded or folded-even site frequency spectra.\n");
    exit(-1);
  }
  l = 0.;
  x = sfs->numPol + sfs->nullCount;
  for(int r = 1; r < max; r++) {
    if(sfs->type == UNFOLDED)
      e = expG(ps->N, sfs->u, ps->m, ps->n, x, ps->k, r);
    else
      e = expF(ps->N, sfs->u, ps->m, ps->n, x, ps->k, r);
    if(!sfs->x[r-1])
      l += sfs->f[r-1] * log(e) - e;
  }
  if(sfs->type == UNFOLDED)
    e = expG(ps->N, sfs->u, ps->m, ps->n, x, ps->k, 0);
  else
    e = expF(ps->N, sfs->u, ps->m, ps->n, x, ps->k, 0);
  l += sfs->nullCount * log(e) - e;
  
  return l;
}

int newtonComp(Sfs *sfs, PopSizes *ps, Args *args) {
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;
  int status;
  size_t i, iter = 0;
  Rparams *p = newRparams(sfs, ps);
  gsl_multiroot_function f; 
  gsl_vector *x = gsl_vector_alloc(ps->m);

  if(sfs->type == UNFOLDED)
    f.f = &unfolded;
  else if(sfs->type == FOLDED_EVEN)
    f.f = &folded;
  else {
    fprintf(stderr, "ERROR[epos]: Can only deal with unfolded and folded-even site frequency spectra.\n");
    exit(-1);
  }
  f.n = ps->m;
  f.params = p;
  for(i = 0; i < ps->m; i++) {
    gsl_vector_set(x, i, ps->iniN[i]);
  }
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
      return status;
    }
    status = gsl_multiroot_test_residual(s->f, 1e-7);
  } while (status == GSL_CONTINUE && iter < 1000);
  for(i = 0; i < ps->m; i++) {
    ps->N[i] = gsl_vector_get(s->x, i);
  }
  gsl_multiroot_fsolver_free(s);
  gsl_vector_free(x);
  free(p);

  return status;
}

int newton(Sfs *sfs, PopSizes *ps, Args *args) {
  int status;

  status = newtonComp(sfs, ps, args);

  return status;
}
