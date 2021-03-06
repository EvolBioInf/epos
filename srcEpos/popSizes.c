/***** popSizes.c *********************************
 * Description: Computation of population sizes
 * based on the mathematics layed out in Peter
 * Pfaffelhuber's memo dated November 8, 2018.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Dec 18 11:19:32 2017
 **************************************************/
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf.h>
#include "eprintf.h"
#include "popSizes.h"
#include "util.h"
#include "newton.h"

/* expG returns the expected value of the r-th entry 
 * in the unfolded site frequency spectrum given r > 0
 * as described in equation (S4).
 */
double expGr(PopSizes *ps, Sfs *sfs, int r){
  int n     = sfs->n;
  int m     = ps->m;
  long l    = sfs->l;
  int *k    = ps->k;
  double u  = sfs->u;
  double *N = ps->N;

  double s = 0.;
  for(int i = 1; i <= m; i++) {
    int x = max(n - k[i] + 1, 0);
    double a = binomial(x, r);
    x = max(n - k[i+1] + 1, 0);
    double b = binomial(x, r);
    s += N[i] * (a - b);
  }
  s *= 4. * u * l / (double)r / binomial(n - 1, r);

  return s;
}

/* expG returns the expected value of the r-th entry 
 * in the unfolded site frequency spectrum for all r
 * as specified in equation (S1).
 */
double expG(PopSizes *ps, Sfs *sfs, int r) {
  if(r == 0){
    double s = 0.;
    for(int i = 1; i <= sfs->n - 1; i++)
      s += expGr(ps, sfs, i);
    s = sfs->l - s;
    return s;
  } else
    return expGr(ps, sfs, r);
}

/* expF returns the expected value of the r-th entry 
 * in the folded site frequency spectrum for all r.
 */
double expF(PopSizes *ps, Sfs *sfs, int r) {
  if(r == 0)
    return expG(ps, sfs, r);
  else
    return expG(ps, sfs, r) + expG(ps, sfs, sfs->n - r);
}

void freePopSizes(PopSizes *ps) {
  free(ps->k);
  free(ps->N);
  free(ps);
}

PopSizes *newPopSizes(Sfs *sfs){
  PopSizes *ps = (PopSizes *)emalloc(sizeof(PopSizes));
  ps->m = 1;
  ps->N = (double *)emalloc(sfs->n       * sizeof(double));
  for(int i = 0; i < sfs->n; i++)
    ps->N[i] = 0.;
  ps->k = (int *)   emalloc((sfs->n + 1) * sizeof(int));
  ps->k[1] = 2;
  ps->k[2] = sfs->n + 1;
  ps->l = 0.;

  return ps;
}

/* dSquared computs the goodness of fit measure on p. 442 of 
 * Lapierre et al. (2017). Accuracy of demographic inferences from
 * the site frequency spectrum: The case of the Yoruba population.
 * Genetics, 206, 439-449. In addition, it computes the expected
 * site frequency spectrum.
 */
void dSquared(PopSizes *ps, Sfs *sfs) {
  double e, o;

  for(int r = 1; r <= sfs->a; r++) {
    if(sfs->G[r] < 0)
      continue;
    if(sfs->f)
      e = expF(ps, sfs, r);
    else
      e = expG(ps, sfs, r);
    sfs->E[r] = e;
    sfs->e += e;
    sfs->o += sfs->G[r];
  }
  for(int r = 1; r <= sfs->a; r++) {
    if(sfs->G[r] < 0)
      continue;
    e = sfs->E[r] / sfs->e;
    o = sfs->G[r] / sfs->o;
    double x = e - o;
    sfs->d += x * x / e;
  }
}

double logLik(PopSizes *ps, Sfs *sfs) {
  double e, l = 0.;

  if (sfs->c == 0) {
    for (int i = 0; i <= sfs->a; i++)
      if (sfs->G[i] != -1)
	sfs->c += gsl_sf_lnfact(sfs->G[i]);
  }
  for(int r = 0; r <= sfs->a; r++) {
    if(sfs->G[r] < 0)
      continue;
    if(sfs->f) {
      e = expF(ps, sfs, r);
    } else {
      e = expG(ps, sfs, r);
    }
    l += (double)sfs->G[r] * log(e) - e;
  }
  l -= sfs->c;

  return l;
}
