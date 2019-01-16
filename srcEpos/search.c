/***** search.c ***********************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Fri Nov 23 15:00:29 2018
 **************************************************/
#include <stdlib.h>
#include <math.h>
#include "sfs.h"
#include "popSizes.h"
#include "interface.h"
#include "search.h"
#include "newton.h"
#include "eprintf.h"

int *nextGreedy(int m, int n, int *start, short setup);
int *nextExhaustive(int n, int m, int *k, short setup);

void cpK(int *k1, int *k2, int m) {
  for(int i = 1; i <= m + 1; i++)
    k2[i] = k1[i];
}

int *nextConfig(int m, int n, int *start, Args *args, short setup) {
  if(m > args->E)
    return nextGreedy(m, n, start, setup);
  else
    return nextExhaustive(m, n, start, setup);
}

void printConfig(int *k, int m, double ll) {
  printf("#m = %d; maximum Log(Likelihood): %f\t", m, ll);
  printf("{%d", k[1]);
  for(int i = 2; i <= m; i++)
    printf(", %d", k[i]);
  printf("}\n");
}

/* compPopSizes computes population sizes and returns their log-likelihood 
 * obtained from cross-validation. 
 */
double compPopSizes(int *kd, int m, Sfs *sfs, PopSizes *ps, Args *args, SfsSet *ss) {
  /* set the new configuration of breakpoints, kd */
  for(int i = 1; i <= m ; i++)
    ps->k[i] = kd[i];
  ps->k[m + 1] = sfs->n + 1;
  ps->m = m;
  /* analyze full data set */
  newton(sfs, ps, args);
  if(ss == NULL)
    return ps->l;
  /* obtain likelihood from cross-validation */
  double l = 0.;
  for(int i = 0; i < args->x; i++) {
    newton(ss->train[i], ps, args);
    double x = logLik(ps, ss->test[i]);
    l += x;
  }
  l /= (double)args->x;

  return l;
}

void freeSearch(int *ka, int *k, int *kp, SfsSet *ss) {
  free(ka);
  free(k);
  free(kp);
  freeSfsSet(ss);
}

PopSizes *searchLevels(Sfs *sfs, Args *args, gsl_rng *r) {
  int *kd, *ka, *kp, *k; /* arrays of levels   */
  int m;
  short improved;

  /* initialize search */
  PopSizes *ps = newPopSizes(sfs);
  ka = (int *)emalloc((sfs->n + 1) * sizeof(int));
  k  = (int *)emalloc((sfs->n + 1) * sizeof(int));
  kp = (int *)emalloc((sfs->n + 1) * sizeof(int));
  cpK(ps->k, ka, ps->m);
  cpK(ps->k,  k, ps->m);
  cpK(ps->k, kp, ps->m);
  SfsSet *ss = splitSfs(sfs, args, r);
  double l = compPopSizes(ka, ps->m, sfs, ps, args, ss);
  double la = l;
  printConfig(k, 1, l);
  /* iterate over the possible number of levels, n */
  for(m = 2; m <= sfs->n; m++) {
    kd = nextConfig(m, sfs->n, k, args, 1);
    improved = 0;
    while((kd = nextConfig(m, sfs->n, k, args, 0)) != NULL) {
      double ld  = compPopSizes(kd, m, sfs, ps, args, ss);
      double ldd = ld;
      if(args->x > 1) {
	ldd = compPopSizes(kd, m, sfs, ps, args, NULL); /* ensure the full data set also gives a sensible result */
      }
      if(ld > la && !isnan(ldd)) {
	improved = 1;
	la = ld;
	cpK(kd, ka, m);
      }
    }
    if(la <= l + args->c) { /* no significant improvement, quit search */
      if(improved)
	printConfig(ka, m, la);
      else
	printf("#m = %d; no improvement\n", m);
      compPopSizes(kp, m - 1, sfs, ps, args, NULL);
      freeSearch(ka, k, kp, ss);
      return ps;
    }
    cpK(ka, k, m);
    cpK(ka, kp, m);
    l = la;
    printConfig(k, m, l);
  }
  compPopSizes(kp, m - 1, sfs, ps, args, NULL);
  freeSearch(ka, k, kp, ss);

  return ps;
}
