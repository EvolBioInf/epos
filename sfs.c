/***** sfs.c **************************************
 * Description: Deal with site frequency spectrum.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Oct 23 10:30:27 2017
 **************************************************/
#include <stdlib.h>
#include "sfs.h"
#include "eprintf.h"
#include "tab.h"
#include "sfs.h"

Sfs *newSfs(int n, int type){
  Sfs *sfs;
  int i;

  sfs = (Sfs *)emalloc(sizeof(Sfs));
  sfs->n = n;
  sfs->f = (double *)emalloc(n*sizeof(double));
  for(i=0; i<n; i++)
    sfs->f[i] = 0.;
  sfs->arr = NULL;
  sfs->type = type;
  sfs->nullCount = 0.;
  sfs->numPol = 0.;
  sfs->u = 0.;

  return sfs;
}

void freeSfs(Sfs *sfs){
  if(sfs == NULL)
    return;
  free(sfs->f);
  if(sfs->arr != NULL)
    free(sfs->arr);
  free(sfs);
}

/* getSfs obtains the next SFS from an open file */
Sfs *getSfs(FILE *fp, Args *args){
  char *line;
  Sfs *sfs;
  int i, type;
  char body;
  double j;

  if(args->U){
    type = UNFOLDED;
  }else{
    type = FOLDED_EVEN;
  }
  sfs = newSfs(0, type);
  body = 0;
  while((line = tabGetLine(fp)) != NULL){
    if(tabNfield() >= 2){
      i = atoi(tabField(0));
      j = atof(tabField(1));
    }else{
      i = 0;
      j = 0;
    }
    /* look for comment and blank lines */
    if(line[0] == '#' || tabNfield() < 2){
      if(body)
	break;
      else
	continue;
    }else if(i == 0 && j > 0){
      body = 1;
      sfs->nullCount = j;
    }else{
      body = 1;
      sfs->f = (double *)erealloc(sfs->f, (sfs->n+1)*sizeof(double));
      sfs->f[sfs->n] = (double)atof(tabField(1));
      sfs->numPol += sfs->f[sfs->n];
      sfs->n++;
    }
  }
  tabFree();
  if(!body) {
    freeSfs(sfs);
    return NULL;
  }
  sfs->u = args->u;
  if(sfs->nullCount == 0){
    if(args->l == 0) {
      fprintf(stderr, "ERROR[epos]: Please include either the zero-class in the SFS or enter the sequence length via the -l option.\n");
      exit(-1);
    }
    sfs->nullCount = args->l - sfs->numPol;
  }
  if(args->U){ /* unfolded */
    sfs->n++;
    sfs->type = UNFOLDED;
  }else{        /* folded/even */
    sfs->n *= 2;
    sfs->type = FOLDED_EVEN;
  }

  return sfs;
}

void iniBoot(Sfs *sfs) {
  int i, j, x, nc, n;

  sfs->arr = (int *)emalloc((sfs->numPol + sfs->nullCount)*sizeof(int));
  n = sfs->n;
  if(sfs->type == FOLDED_EVEN)
    n /= 2;
  nc = sfs->nullCount;
  for(i=0;i<nc;i++)
    sfs->arr[i] = 0;
  x = nc;
  for(i=0;i<n-1;i++)
    for(j=0;j<sfs->f[i];j++)
      sfs->arr[x++] = i+1;
  if(x != sfs->numPol + sfs->nullCount){
    printf("ERROR in bootstrapping; sfs->numPol: %d; sfs->nullCount: %d; x: %d; expected(x): %d. Is this an unfolded SFS? If so, please use -U.\n", (int)sfs->numPol, (int)sfs->nullCount, x, (int)(sfs->numPol+sfs->nullCount));
    exit(-1);
  }
}

Sfs *bootstrapSfs(Sfs *sfs, gsl_rng *rand, Args *args){
  int i, j, x, y, type;
  Sfs *bootSfs;
  
  if(args->U)
    type = UNFOLDED;
  else
    type = FOLDED_EVEN;
  bootSfs = newSfs(sfs->n, type);
  bootSfs->u = args->u;
  x = sfs->numPol + sfs->nullCount;
  for(i=0; i<x; i++){
    y = gsl_rng_uniform(rand) * x;
    j = sfs->arr[y];
    if(j){
      bootSfs->f[j-1]++;
      bootSfs->numPol++;
    }else
      bootSfs->nullCount++;
  }

  return bootSfs;
}
