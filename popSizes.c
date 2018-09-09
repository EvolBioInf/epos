/***** popSizes.c *********************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Dec 18 11:19:32 2017
 **************************************************/
#include <stdlib.h>
#include <omp.h>
#include "eprintf.h"
#include "popSizes.h"
#include "printMatrix.h"
#include "util.h"
#include "newton.h"

void findIniN(PopSizes *ps, Sfs *sfs) {
  int id;
  char found;

  id = -1;
  if(ps->m == 1) {
    ps->iniN[0] = watterson(sfs);
    return;
  }
  for(int i = 1; i < ps->m; i++) {
    found = 0;
    for(int j = 1; j < ps->m - 1; j++) {
      if(ps->k[i] == ps->prevK[j]) {
	found = 1;
      }
    }
    if(!found) {
      id = i;
      break;
    }
  }
  for(int i = 0; i < id; i++)
    ps->iniN[i] = ps->N[i];
  for(int i = id; i < ps->m; i++)
    ps->iniN[i] = ps->N[i - 1];
}

int negPopSizes(PopSizes *ps){
  int i;

  for(i=0; i<ps->m; i++)
    if(ps->N[i] < 1)
      return 1;

  return 0;
}

void freePopSizes(PopSizes *ps){
  free(ps->k);
  free(ps->N);
  free(ps->prevK);
  free(ps->iniN);
  free(ps);
}

PopSizes *newPopSizes(Sfs *sfs){
  PopSizes *ps;
  int i, n;

  ps = (PopSizes *)emalloc(sizeof(PopSizes));
  ps->N = (double *)emalloc(sfs->n * sizeof(double));
  ps->iniN = (double *)emalloc(sfs->n * sizeof(double));
  ps->k = (int *)emalloc((sfs->n+1) * sizeof(int));
  ps->prevK = (int *)emalloc((sfs->n+1) * sizeof(int));
  n = sfs->n;
  for(i=0; i<n; i++){
    ps->N[i]    = 0;
    ps->k[i]    = 0;
    ps->iniN[i] = 0;
  }
  ps->k[0] = n+1;
  ps->prevK[0] = n+1;
  ps->m = 0;
  ps->prevM = 0;

  return ps;
}

int cmpInt(const void * a, const void * b){
  return *(int*)a - *(int*)b;
}

void addTestK(PopSizes *ps, int k){
  int i;

  for(i=0;i<ps->prevM+1;i++)
    ps->k[i] = ps->prevK[i];
  ps->k[ps->prevM+1] = k;
  ps->m = ps->prevM + 1;
  qsort(ps->k,ps->m+1,sizeof(int),cmpInt);
}

void addK(PopSizes *ps, int k){
  int i;

  addTestK(ps,k);
  for(i=0;i<ps->m+1;i++)
    ps->prevK[i] = ps->k[i];
  ps->prevM++;
}

void restoreK(PopSizes *ps){
  int i;

  for(i=0;i<ps->prevM+1;i++)
    ps->k[i] = ps->prevK[i];
  ps->m--;
}

double psi(PopSizes *ps, Sfs *sfs){
  return logLik(ps, sfs);
}

int compPopSizes(Sfs *sfs, PopSizes *ps) {
  return newton(sfs, ps);
}

double testK(Sfs *sfs, PopSizes *ps, int k){
  int status;

  addTestK(ps, k);
  if((status = compPopSizes(sfs, ps)) > 0)
    return DBL_MIN;
  if(negPopSizes(ps))
    return DBL_MIN;

  return logLik(ps, sfs);
}

PopSizes *copyPopSizes(PopSizes *ps){
  PopSizes *cps;
  int i;

  cps = (PopSizes *)emalloc(sizeof(PopSizes));
  cps->N = (double *)emalloc((ps->m+2) * sizeof(double));
  cps->k = (int *)emalloc((ps->m+2) * sizeof(int));
  cps->prevK = (int *)emalloc((ps->m+2) * sizeof(int));
  cps->n = ps->n;
  cps->m = ps->m;
  cps->prevM = ps->prevM;
  cps->psi = ps->psi;

  for(i=0; i<cps->m+1; i++){
    cps->N[i] = ps->N[i];
    cps->k[i] = ps->k[i];
    cps->prevK[i] = ps->prevK[i];
  }

  return cps;
}

int getNextLevel(Sfs *sfs, PopSizes *ps, int *avail){
  double p, currMinPsi;
  int i, minK;

  currMinPsi = DBL_MIN;
  minK = 0;
  for(i=3; i<=sfs->n; i++){
    if(avail[i]){
      p = testK(sfs, ps, i);
      if(p > currMinPsi){
	minK = i;
	currMinPsi = p;
      }
    }
  }
  avail[minK] = 0;

  return minK;
}

PopSizes *getPopSizes(Sfs *sfs){
  int i, l, n, minK, *avail;
  double prevMinPsi, currMinPsi, change;
  PopSizes *ps;

  ps = newPopSizes(sfs);
  n = sfs->n;
  ps->n = sfs->n;
  /* add first entry, single pop size for entire coalescent */
  minK = 2;
  testK(sfs, ps, minK);
  addK(ps, minK);
  compPopSizes(sfs, ps);
  ps->psi = psi(ps, sfs);
  prevMinPsi = ps->psi;
  /* find remaining entries */
  /* keep track of which entries are still available */
  avail = (int *)emalloc((n+2)*sizeof(int));
  for(i=3; i<=n; i++)
    avail[i] = 1;
  /* iterate over the remaining number of possible population sizes */
  for(i=3; i<=n; i++){
    l = getNextLevel(sfs, ps, avail);
    if(l)
      currMinPsi = testK(sfs, ps, l);
    else
      currMinPsi = DBL_MIN;
    change = currMinPsi - prevMinPsi;
    if(change < 2){
      restoreK(ps);
      compPopSizes(sfs, ps);
      ps->psi = psi(ps, sfs);
      break;
    }else{
      addK(ps, l);
      compPopSizes(sfs, ps);
      currMinPsi = psi(ps, sfs);
    }
    prevMinPsi = currMinPsi;
  }
  free(avail);

  return ps;
}
