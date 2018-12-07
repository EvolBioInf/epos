/***** ages.c *************************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Thu Nov 29 12:31:24 2018
 **************************************************/
#include <stdlib.h>
#include "popSizes.h"
#include "eprintf.h"
#include "ages.h"
#include "util.h"

void printAges(Ages *a) {
  printf("#r\tA[r]\tV(A[r])\tP[r]\n");	
  for(int r = 1; r < a->n; r++)
    printf("%d\t%g\t%g\t%g\n", r, a->a[r], a->v[r], a->s[r]);
}

/* numAaa computes the numerator in the equation for the average age of an allele */
double numAaa(Ages *a, int r) {
  double *N = a->N;
  int     n = a->n;
  double nu = 0.;

  for(int k = 2; k <= n; k++) {
    double s = 0.;
    for(int l = k; l <= n; l++)
      s += N[l] / (double) l / (double) (l-1);
    nu += N[k] * binomial(n-k, r-1) * s;
  }
  nu *= 4;

  return nu;
}

/* aaa computes the average age of an allele */
double aaa(Ages *a, int r) {
  double *N = a->N;
  int     n = a->n;
  /* Numerator */
  double nu = numAaa(a, r);
  /* Denomiator */
  double de = 0.;
  for(int k = 2; k <= n; k++)
    de += N[k] * binomial(n-k, r-1);

  return nu / de;
}

/* varAaa computes the variance of the age of an allele
 * Reference: Peter's memo of Dec. 3, 2018
 */
double varAaa(Ages *a, int r) {
  double s1, s2;
  double *N = a->N;
  int     n = a->n;

  double nu = 0.;                 /* numerator */
  for(int k = 2; k <= n; k++) {
    s1 = 0.;
    for(int l = k; l <= n; l++) {
      s2 = 0.;
      for(int m = l; m <= n; m++) {
	s2 += N[m] / (double) m / (double) (m - 1);
      }
      s1 += N[l] / (double) l / (double) (l - 1) * s2;
    }
    nu += N[k] * binomial(n - k, r - 1) * s1 * s2;
  }
  nu *= 32.;

  double de = 0.;                 /* denominator */
  for(int k = 2; k <= n; k++)
    de += N[k] * binomial(n - k, r - 1);
  double x2 = nu / de;

  double x = aaa(a, r);
  double v = x2 - x * x;

  return v;
}

/* asa computes the average size of an allele */
double asa(Ages *a, int r) {
  double *N = a->N;
  int     n = a->n;
  /* Numerator */
  double nu = 0.;
  for(int k = 2; k <= n; k++) {
    double s = 0.;
    for(int l = k; l <= n; l++)
      s += 4. * N[l] * N[l] / (double) l / (double) (l-1);
    nu += N[k] * binomial(n-k, r-1) * s;
  }
  /* Denomiator */
  double de = numAaa(a, r);

  return nu / de;
}

void freeAges(Ages *a) {
  free(a->N);
  free(a->a);
  free(a->v);
  free(a->s);
  free(a);
}

Ages *newAges(PopSizes *ps) {
  Ages *a = (Ages *)emalloc(sizeof(Ages));
  int n = ps->k[ps->m + 1] - 1;
  a->n = n;
  a->N = (double *)emalloc((n + 1) * sizeof(double));
  a->a = (double *)emalloc((n + 1) * sizeof(double));
  a->v = (double *)emalloc((n + 1) * sizeof(double));
  a->s = (double *)emalloc((n + 1) * sizeof(double));

  for(int i = 1; i <= ps->m; i++)
    for(int j = ps->k[i]; j < ps->k[i+1]; j++)
      a->N[j] = ps->N[i];

  return a;
}

Ages *compAges(PopSizes *ps) {
  Ages *a = newAges(ps);
  for(int r = 1; r < a->n; r++) {
    a->a[r] = aaa(a, r);
    a->v[r] = varAaa(a, r);
    a->s[r] = asa(a, r);
  }
  return a;
}
