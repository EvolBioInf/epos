/***** util.c *************************************
 * Description: Equations taken from Peter Pfaffel-
 *   huber's write-up dated December 19, 2017.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Fri Dec 15 10:36:55 2017
 **************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_sf_gamma.h>
#include <math.h>
#include "util.h"
#include "eprintf.h"

double **bin = NULL;
int nn;

void iniBinom(int n) {
  double x, y;

  nn = n;

  bin = (double  **)emalloc((n + 1) * sizeof(double *));
  for(int i = 0; i <= n; i++) {
    bin[i] = (double *)emalloc((n + 1) * sizeof(double));
    for(int j = 0; j <= n; j++) {
      if(i >= j) {
	x = gsl_sf_lnfact(i) - (gsl_sf_lnfact(j) + gsl_sf_lnfact(i-j));
	y = exp(x);
	y = round(y);
	bin[i][j] = y;
      }
    }
  }
}

void freeBinom() {
  for(int i = 0; i <= nn; i++)
    free(bin[i]);
  free(bin);
  bin = NULL;
}

double binomial(int n, int k){

  if(n < k)
    return 0;
  else
    return bin[n][k];
  
}

/* fi: Equation (5) */
double fi(int i, int n, int r, int *k){
  double a, b;

  a = binomial(n-k[i-1]+1, r);
  b = binomial(n-k[i]  +1, r);

  return a - b;
}

/* gi: Equation in the middle of p. 5 */
double gi(int i, int n, int r, int *k){
  double a, b;

  a = fi(i, n, r, k);
  b = fi(i, n, n-r, k);

  return a + b;
}

/* hi: Equation at the bottom of p. 6 */
double hi(int i, int n, int r, int *k){
  if(r < n/2)
    return gi(i, n, r, k);
  else if(r == n/2)
    return fi(i, n, r, k);
  else{
    fprintf(stderr, "ERROR [util.hi]: r > n/2!\n");
    exit(-1);
  }
}

void printSfsStats(Sfs *sfs){
  printf("#Polymorphic sites surveyed:\t%d\n", sfs->p);
  printf("#Monomorphic sites surveyed:\t%d\n", sfs->G[0]);
}

void printTimes(PopSizes *ps, Sfs *sfs){
  int i, k;
  double t;

  k = sfs->n;
  if(ps->l > 0)
    printf("#Final Log(Likelihood):          %f\n", ps->l);
  else
    printf("#Final Log(Likelihood):          %f\n", 0.0);
  printf("#Level\tT[Level]\tN[Level]\n");
  t = 0.;
  for(i=ps->m;i>=1;i--){
    for(;k>=ps->k[i];k--)
      t += 4. / (double)k / (double)(k-1) * ps->N[i];
    printf("%d\t%.2e\t%.2e\n", k+1, t, ps->N[i]);
    if(k + 1 == 2)
      break;
  }
}

/* watterson: Using Watterson's estimator of N */
double watterson(Sfs *sfs){
  int i;
  double s, l, w;

  s = 0.;
  for(i=1; i<sfs->n; i++)
    s += 1./i;
  l = sfs->G[0] + sfs->p;
  w = sfs->p / s / 4. / sfs->u / l;

  return w;
}

int max(int a, int b) {
  if(a > b)
    return a;
  else
    return b;
}
