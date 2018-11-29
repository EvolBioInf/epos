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
