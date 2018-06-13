/***** popSizes.c *********************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Dec 18 11:19:32 2017
 **************************************************/
#include <stdlib.h>
#include <omp.h>
#include "eprintf.h"
#include "popSizes.h"
#include "unfolded.h"
#include "foldedE.h"
#include "printMatrix.h"
#include "util.h"

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
  free(ps);
}

PopSizes *newPopSizes(Sfs *sfs){
  PopSizes *ps;
  int i, n;

  ps = (PopSizes *)emalloc(sizeof(PopSizes));
  ps->N = (double *)emalloc(sfs->n * sizeof(double));
  ps->k = (int *)emalloc((sfs->n+1) * sizeof(int));
  ps->prevK = (int *)emalloc((sfs->n+1) * sizeof(int));
  n = sfs->n;
  for(i=0; i<n; i++){
    ps->N[i] = 0.;
    ps->k[i] = 0;
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
  if(sfs->type == UNFOLDED){
    return unfoldedPsi(ps, sfs);
  }
  else if(sfs->type == FOLDED_EVEN)
    return foldedEpsi(ps, sfs);
  else if(sfs->type == FOLDED_ODD){
    fprintf(stderr, "ERROR[psi]: FoldedOdd not implemented yet.\n");
    exit(-1);
  }else{
    fprintf(stderr, "ERROR[psi]: Requesting SFS-type other than unfolded|foldedOdd|foldedEven.\n");
    exit(-1);
  }
}

gsl_vector *getResVect(Sfs *sfs, PopSizes *ps){
  if(sfs->type == UNFOLDED)
    return unfoldedGetResVect(sfs, ps);
  else if(sfs->type == FOLDED_EVEN)
    return foldedEgetResVect(sfs, ps);
  else if(sfs->type == FOLDED_ODD){
    fprintf(stderr, "ERROR[getResVect]: FoldedOdd not implemented yet.\n");
    exit(-1);
  }else{
    fprintf(stderr, "ERROR[getResVect]: Requesting SFS-type other than unfolded|foldedOdd|foldedEven.\n");
    exit(-1);
  }
}

gsl_matrix *getCoeffMat(Sfs *sfs, PopSizes *ps, Args *args){
  if(sfs->type == UNFOLDED)
    return unfoldedGetCoeffMat(sfs, ps, args);
  else if(sfs->type == FOLDED_EVEN)
    return foldedEgetCoeffMat(sfs, ps, args);
  else if(sfs->type == FOLDED_ODD){
    fprintf(stderr, "ERROR[getCoeffMat]: FoldedOdd not implemented yet.\n");
    exit(-1);
  }else{
    fprintf(stderr, "ERROR[getCoeffMat]: Requesting SFS-type other than unfolded|foldedOdd|foldedEven.\n");
    exit(-1);
  }
}

void freeGsl(gsl_matrix *A, gsl_matrix *LU, gsl_vector *b, gsl_permutation *p, gsl_vector *x, gsl_vector *residual){
  gsl_matrix_free(A);
  gsl_matrix_free(LU);
  gsl_vector_free(b);
  gsl_permutation_free(p);
  gsl_vector_free(x);
  gsl_vector_free(residual);
}

int compPopSizes(Sfs *sfs, PopSizes *ps, Args *args){
  int i, s, status;
  gsl_matrix *LU, *A;
  gsl_vector *b;


  gsl_vector *residual = gsl_vector_alloc(ps->m);
  gsl_vector *x = gsl_vector_alloc(ps->m);
  gsl_permutation *p = gsl_permutation_alloc(ps->m);

  A = gsl_matrix_alloc(ps->m, ps->m);
  b = getResVect(sfs, ps);
  LU = getCoeffMat(sfs, ps, args);
  if(args->p)
    printMatrix(LU, b);
  gsl_matrix_memcpy(A, LU);
  gsl_set_error_handler_off();
  status = gsl_linalg_LU_decomp(LU, p, &s);
  if(status){
    if(args->V)
      printf("#Error1 in linear algebra; returning maximal value for Psi\n");
    freeGsl(A, LU, b, p, x, residual);
    return status;
  }
  status = gsl_linalg_LU_solve(LU, p, b, x);
  if(status){
    if(args->V)
      printf("#Error2 in linear algebra; returning maximal value for Psi.\n");
    freeGsl(A, LU, b, p, x, residual);
    return status;
  }
  status = gsl_linalg_LU_refine(A, LU, p, b, x, residual);
  if(status){
    if(args->V)
      printf("#Error3 in linear algebra; returning maximal value for Psi\n");
    freeGsl(A, LU, b, p, x, residual);
    return status;
  }
  for(i=0; i<ps->m; i++)
    ps->N[i] = x->data[i];
  /* gsl_set_error_handler(NULL); */
  freeGsl(A, LU, b, p, x, residual);
  return 0;
}

double testK(Sfs *sfs, PopSizes *ps, Args *args, int k){
  addTestK(ps, k);
  if(compPopSizes(sfs, ps, args))
    return DBL_MAX;
  if(negPopSizes(ps) && !args->n){
    if(args->V)
      fprintf(stderr,"#Negative population size\n");
    return DBL_MAX;
  }
  return psi(ps, sfs);
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

int getNextLevel(Sfs *sfs, PopSizes *ps, Args *args, int *avail){
  double p, currMinPsi;
  int i, minK;

  currMinPsi = DBL_MAX;
  minK = 0;
  for(i=3; i<=sfs->n; i++){
    if(avail[i]){
      p = testK(sfs, ps, args, i);
      if(args->V)
	printTimes(ps, sfs);
      if(p < currMinPsi){
	minK = i;
	currMinPsi = p;
      }
    }
  }
  avail[minK] = 0;

  return minK;
}

PopSizes *getPopSizes(Sfs *sfs, Args *args){
  int i, l, n, minK, *avail;
  double prevMinPsi, currMinPsi, change;
  PopSizes *ps;

  ps = newPopSizes(sfs);
  n = sfs->n;
  ps->n = sfs->n;
  /* add first entry, single pop size for entire coalescent */
  minK = 2;
  if(!args->m)
    args->m = n-1;
  testK(sfs, ps, args, minK);
  addK(ps, minK);
  compPopSizes(sfs, ps, args);
  ps->psi = psi(ps, sfs);
  prevMinPsi = ps->psi;
  if(args->V){
    printf("#Watterson: %.2e\n", watterson(sfs));
    printTimes(ps, sfs);
  }
  /* find remaining entries */
  /* keep track of which entries are still available */
  avail = (int *)emalloc((n+2)*sizeof(int));
  for(i=3; i<=n; i++)
    avail[i] = 1;
  /* iterate over the remaining number of possible population sizes */
  for(i=3; i<=n; i++){
    if(ps->m >= args->m){
      ps->psi = prevMinPsi;
      free(avail);
      return ps;
    }
    l = getNextLevel(sfs, ps, args, avail);
    if(l)
      currMinPsi = testK(sfs, ps, args, l);
    else
      currMinPsi = DBL_MAX;
    change = prevMinPsi - currMinPsi;
    if(change <= args->d){
      if(args->V){
	ps->psi = psi(ps, sfs);
	printTimes(ps, sfs);
	printf("#Reverting to previous configuration\n");
      }
      restoreK(ps);
      compPopSizes(sfs, ps, args);
      ps->psi = psi(ps, sfs);
      break;
    }else{
      addK(ps, l);
      compPopSizes(sfs, ps, args);
      currMinPsi = psi(ps, sfs);
    }
    prevMinPsi = currMinPsi;
  }
  free(avail);

  return ps;
}
