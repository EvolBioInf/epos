/***** aaa.c **************************************
 * Description: Average age of an allele. For de-
 *   rivation, see Peter Pfaffelhuber's memo of
 *   August 20, 2018.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Tue Sep 18 08:09:08 2018
 **************************************************/
#include "popSizes.h"
#include "eprintf.h"
#include "util.h"

/* numAaa computes the numerator in the equation for the average age of an allele */
double numAaa(double *N, int n, int r) {
  double nu = 0.;
  for(int k = 2; k <= n; k++) {
    double s = 0.;
    for(int l = k; l <= n; l++) {
      s += 4. * N[l-1] / (double) l / (double) (l-1);
    }
    nu += N[k-1] * binomial(n-k, r-1) * s;
  }
  return nu;
}

/* aaa computes the average age of an allele */
double aaa(double *N, int n, int r) {

  /* Numerator */
  double nu = numAaa(N, n, r);
  /* Denomiator */
  double de = 0.;
  for(int k = 2; k <= n; k++) {
    de += N[k-1] * binomial(n-k, r-1);
  }

  return nu / de;
}

/* asa computes the average size of an allele */
double asa(double *N, int n, int r) {
  /* Numerator */
  double nu = 0.;
  for(int k = 2; k <= n; k++) {
    double s = 0.;
    for(int l = k; l <= n; l++) {
      s += 4. * N[l-1] * N[l-1] / (double) l / (double) (l-1);
    }
    nu += N[k-1] * binomial(n-k, r-1) * s;
  }
  /* Denomiator */
  double de = numAaa(N, n, r);

  return nu / de;
}

void allN(PopSizes *ps) {
  for(int i = 0; i < ps->m; i++) {
    for(int j = ps->k[i]; j <= ps->k[i+1]; j++)
      ps->allN[j-1] = ps->N[i];
  }
}


void compAaa(PopSizes *ps, Sfs *sfs) {
  int n = ps->n;
  int max;

  if(sfs->type == UNFOLDED)
    max = n - 1;
  else
    max = n / 2;
  allN(ps);
  for(int r = 1; r <= max; r++) {
    ps->aaa[r-1] = aaa(ps->allN, n, r);
    ps->asa[r-1] = asa(ps->allN, n, r);
  }
}
