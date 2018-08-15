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

double binomial(int n, int k){
  double x, y;
  
  if(n==k)
    return 1;
  if(n < k)
    return 0;

  x = gsl_sf_lnfact(n) - (gsl_sf_lnfact(k) + gsl_sf_lnfact(n-k));
  y = exp(x);
  y = round(y);
  
  
return y;
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

/* delta: As explained by Peter */
int delta(int x, int y){
  return x==y ? 1 : 0;
}

void printSfsStats(Sfs *sfs){
  printf("#Polymorphic sites surveyed:\t%g\n", sfs->numPol);
  printf("#Monomorphic sites surveyed:\t%g\n", sfs->nullCount);
}

void printTimes(PopSizes *ps, Sfs *sfs){
  int i, k;
  double t;

  k = ps->n;
  printf("#lambda:\t%e\n", sfs->l);
  printf("#Psi:\t\t%e\n", ps->psi);
  printf("#Level\tT[Level]\tN[T]\n");
  t = 0.;
  for(i=ps->m-1;i>=0;i--){
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
  double s;

  s = 0.;
  for(i=1; i<sfs->n; i++)
    s += 1./i;

  return sfs->numPol / s / 4. / sfs->u;
}
