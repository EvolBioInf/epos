/***** popSizes.c *********************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Thu Nov 29 12:22:26 2018
 **************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "tab.h"
#include "eprintf.h"
#include "popSizes.h"

PopSizes *newPopSizes() {
  PopSizes *p = (PopSizes *)emalloc(sizeof(PopSizes));
  p->m = 0;
  p->N = (double *)emalloc(sizeof(double));
  p->k = (int *)emalloc(sizeof(int));
  return p;
}

void freePopSizes(PopSizes *p) {
  free(p->N);
  free(p->k);
  free(p);
}

void reverse(PopSizes *p) {
  int *k    = (int *)    emalloc((p->m + 2) * sizeof(int));
  double *N = (double *) emalloc((p->m + 2) * sizeof(double));
  int j = 1;
  for(int i = p->m; i >= 1; i--) {
    k[j] = p->k[i];
    N[j] = p->N[i];
    j++;
  }
  for(int i = 1; i <= p->m; i++) {
    p->k[i] = k[i];
    p->N[i] = N[i];
  }
  free(N);
  free(k);
}

PopSizes *nextPs(FILE *fp) {
  static short open = 0;
  static short last = 0;
  char *line;

  if(last)
    return NULL;
  PopSizes *ps = newPopSizes();
  while((line = tabGetLine(fp)) != NULL) {
    if(line[0] == '#') {
      if(open) {
	open = 0;
	reverse(ps);
	return ps;
      }
      continue;
    } else {
      if(!open)
	open = 1;
      ps->m++;
      int l    = atoi(tabField(0));
      double s = atof(tabField(2));
      ps->N = (double *)erealloc(ps->N, (ps->m + 2) * sizeof(double));
      ps->k = (int *)   erealloc(ps->k, (ps->m + 2) * sizeof(int));
      ps->k[ps->m] = l;
      ps->N[ps->m] = s;
    }
  }
  last = 1;
  reverse(ps);
  return ps;
}
