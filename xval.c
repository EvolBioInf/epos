/***** xval.c *************************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Thu Mar  8 17:45:00 2018
 **************************************************/
#include <stdio.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_errno.h>
#include "sfs.h"
#include "popSizes.h"
#include "util.h"

void summarizeSfs(Sfs *recipient, Sfs **sfsArr, int n, int excluded){
  int i;

  resetSfs(recipient);
  for(i=0;i<n;i++)
    if(i != excluded)
      addSfs(recipient, sfsArr[i]);

}

double testError(Sfs *trainSfs, Sfs **sfsArr, Args *args, double testLambda){
  int i;
  PopSizes *ps;
  double e, testPsi;
  Sfs *testSfs;

  e = 0;
  for(i=0; i<args->c; i++){
    testSfs = sfsArr[i];
    testSfs->l = testLambda;
    summarizeSfs(trainSfs, sfsArr, args->c, i);
    ps = getPopSizes(trainSfs, args);
    testPsi = psi(ps, testSfs);
    e += testPsi;
    
    freePopSizes(ps);
  }

  return e;
}

/* xvalLambda: Find best lambda through cross-validation. */
double xvalLambda(Sfs *trainSfs, Sfs **sfsArr, Args *args){
  double e, minErr, minL;
  
  minL = 0.;
  trainSfs->l = minL;
  minErr = testError(trainSfs, sfsArr, args, 0.);
  if(args->V)
    printf("#Lambda(m=%d): %e; Error: %e\n", args->m, minL, minErr);
  trainSfs->l += DBL_EPSILON;
  while(trainSfs->l < trainSfs->u){
    e = testError(trainSfs, sfsArr, args, 0.);
    if(args->V)
      printf("#Lambda(m=%d): %e; Error: %e\n", args->m, trainSfs->l, e);
    if(e < minErr){
      minErr = e;
      minL = trainSfs->l;
    }
    trainSfs->l *= 2.;
  }

  return minL;
}

/* xvalM: Find optimal number of levels, m, by cross-validation. */
void xvalM(Sfs *sfs, Args *args, gsl_rng *r){
  int i;
  double currTestError, prevTestError, currLambda, prevLambda;
  Sfs *trainSfs, **sfsArr;

  sfsArr = splitSfs(sfs, r, args);       /* split input sfs */
  trainSfs = newSfs(sfs->n, sfs->type);  /* create empty training SFS */
  trainSfs->l = sfs->l;
  prevLambda = currLambda = sfs->l;
  trainSfs->u = sfs->u;
  prevTestError = DBL_MAX;
  for(i=1; i<=sfs->n; i++){
    args->m = i;
    if(args->L && i>1)
      currLambda = xvalLambda(trainSfs, sfsArr, args);
    trainSfs->l = currLambda;
    currTestError = testError(trainSfs, sfsArr, args, currLambda);
    if(args->V)
      printf("#Avg. cross-validation error, m=%d: %e\n", i, currTestError/args->c);
    if(currTestError >= prevTestError){
      args->m--;
      if(args->L)
	sfs->l = prevLambda;
      break;
    }
    prevTestError = currTestError;
    prevLambda = currLambda;
  }
  for(i=0; i<args->c; i++)
    freeSfs(sfsArr[i]);
  free(sfsArr);
  freeSfs(trainSfs);
}
