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

double aaa(double *N, int n, int r) {
  double de = 0.;

  /* Denominator */
  for(int k = 2; k <= n; k++) {
    double s = 0.;
    for(int l = k; l <= n; l++) {
      s += 4. * N[l-1] / (double) l / (double) (l-1);
    }
    de += N[k-1] * binomial(n-k, r-1) * s;
  }

  /* Numerator */
  double nu = 0.;
  for(int k = 2; k <= n; k++) {
    nu += N[k-1] * binomial(n-k, r-1);
  }

  return de / nu;
}

void compAaa(PopSizes *ps, Sfs *sfs) {
  int n = ps->n;
  int max;

  if(sfs->type == UNFOLDED)
    max = n - 1;
  else
    max = n / 2;

  for(int r = 1; r <= max; r++) {
    ps->aaa[r-1] = aaa(ps->N, n, r);
  }
}
