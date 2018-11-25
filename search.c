/***** search.c ***********************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Fri Nov 23 15:00:29 2018
 **************************************************/
#include <stdlib.h>
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

void printConfig(int *k, int m) {
  printf("%d", k[1]);
  for(int i = 2; i <= m + 1; i++)
    printf(" %d", k[i]);
  printf("\n");
}

/* compPopSizes computes population sizes and returns their log-likelihood */
double compPopSizes(int *kd, int m, Sfs *sfs, PopSizes *ps, Args *args) {
  for(int i = 1; i <= m ; i++)
    ps->k[i] = kd[i];
  ps->k[m + 1] = sfs->n + 1;
  ps->m = m;
  newton(sfs, ps, args);
  return ps->l;
}

PopSizes *searchLevels(Sfs *sfs, Args *args) {
  int *kd, *ka, *k; /* arrays of levels   */
  int m;

  /* initialize search */
  PopSizes *ps = newPopSizes(sfs);
  ka = (int *)emalloc((sfs->n + 1) * sizeof(int));
  k  = (int *)emalloc((sfs->n + 1) * sizeof(int));
  newton(sfs, ps, args);
  double l = ps->l;
  double la = l;
  cpK(ps->k, ka, ps->m);
  cpK(ps->k,  k, ps->m);
  /* iterate over the possible number of levels, n */
  for(m = 2; m <= sfs->n; m++) {
    kd = nextConfig(m, sfs->n, ka, args, 1);
    while((kd = nextConfig(m, sfs->n, ka, args, 0)) != NULL) {
      double ld = compPopSizes(kd, m, sfs, ps, args);
      if(ld > la) {
	la = ld;
	cpK(kd, ka, m);
      }
    }
    if(la <= l + 2) {
      compPopSizes(k, m - 1, sfs, ps, args);
      return ps;
    }
    cpK(ka, k, m);
    l = la;
  }
  compPopSizes(k, m - 1, sfs, ps, args);
  free(k);
  free(ka);

  return ps;
}
