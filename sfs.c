/***** sfs.c **************************************
 * Description: Deal with site frequency spectrum.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Oct 23 10:30:27 2017
 **************************************************/
#include <stdlib.h>
#include "sfs.h"
#include "eprintf.h"
#include "tab.h"
#include "gsl_rng.h"
#include "sfs.h"

Sfs *newSfs(int n, int type){
  Sfs *sfs;
  int i;

  sfs = (Sfs *)emalloc(sizeof(Sfs));
  sfs->n = n;
  sfs->f = (double *)emalloc(n*sizeof(double));
  for(i=0; i<n; i++)
    sfs->f[i] = 0;
  sfs->type = type;
  if(type == UNFOLDED)
    sfs->isFolded = 0;
  else
    sfs->isFolded = 1;
  sfs->nullCount = 0.;
  sfs->numPol = 0.;
  sfs->u = 0.;
  sfs->l = 0.;

  return sfs;
}

Sfs *bootstrapSfs(Sfs *sfs, gsl_rng *rand, Args *args){
  int *arr;
  int i, j, n, nc, x, y, type;
  Sfs *bootSfs;
  
  if(args->U)
    type = UNFOLDED;
  else
    type = FOLDED_EVEN;
  bootSfs = newSfs(sfs->n, type);
  bootSfs->u = args->u * (sfs->numPol + sfs->nullCount);
  arr = (int *)emalloc((sfs->numPol + sfs->nullCount)*sizeof(int));
  n = sfs->n;
  if(sfs->isFolded)
    n /= 2;
  x = 0;
  for(i=0;i<n;i++)
    for(j=0;j<sfs->f[i];j++)
      x++;
  nc = sfs->nullCount;
  for(i=0;i<nc;i++)
    arr[i] = 0;
  x = nc;
  for(i=0;i<n;i++)
    for(j=0;j<sfs->f[i];j++)
      arr[x++] = i+1;
  if(x != sfs->numPol + sfs->nullCount){
    printf("ERROR in bootstrapping; sfs->numPol: %d; sfs->nullCount: %d; x: %d; expected(x): %d\n", (int)sfs->numPol, (int)sfs->nullCount, x, (int)(sfs->numPol+sfs->nullCount));
    exit(-1);
  }
  for(i=0; i<x; i++){
    y = gsl_rng_uniform(rand) * x;
    j = arr[y];
    if(j){
      bootSfs->f[j-1]++;
      bootSfs->numPol++;
    }else
      bootSfs->nullCount++;
  }
  free(arr);

  return bootSfs;
}

/* shuffleArr: Shuffle the entries in array arr. */
void shuffleArr(int *arr, int n, gsl_rng *rand){
  int r, tmp, i;
    
  for(i=n-1; i>0; i--){
    r = (int)(gsl_rng_uniform(rand) * (i+1));
    tmp = arr[r];
    arr[r] = arr[i];
    arr[i] = tmp;
  }
}

void freeSfs(Sfs *sfs){
  if(sfs == NULL)
    return;
  free(sfs->f);
  free(sfs);
}

Sfs *getSfs(FILE *fp, Args *args){
  char *line;
  Sfs *sfs;
  int i, type;
  double j;

  if(args->U)
    type = UNFOLDED;
  else
    type = FOLDED_EVEN;
  sfs = newSfs(0, type);
  while((line = tabGetLine(fp)) != NULL){
    if(tabNfield() >= 2){
      i = atoi(tabField(0));
      j = atof(tabField(1));
    }else{
      i = 0;
      j = 0;
    }
    if(line[0] == '#'){
      continue;
    }else if(i == 0 && j > 0){
      sfs->nullCount = j;
    }else{
      sfs->f = (double *)erealloc(sfs->f, (sfs->n+1)*sizeof(double));
      sfs->f[sfs->n] = (double)atof(tabField(1));
      sfs->numPol += sfs->f[sfs->n];
      sfs->n++;
    }
  }
  j = 0;
  for(i=0; i<sfs->n; i++)
    j += sfs->f[i];
  tabFree();
  sfs->u = args->u;
  /* if monomorphic sites are included in the data, compute sample-wide mutation rate */
  if(sfs->nullCount > 0)
    sfs->u *= (sfs->nullCount + sfs->numPol);
  else
    sfs->u = args->u;
  if(args->U){ /* unfolded */
    sfs->n++;
    sfs->type = UNFOLDED;
  }else{        /* folded/even */
    sfs->n = 2*sfs->n;
    sfs->type = FOLDED_EVEN;
  }
  /* compute lambda as function of mutation rate? */
  if(args->l < 0)
    sfs->l = args->f*sfs->u;
  else
    sfs->l = args->l;

  return sfs;
}

void printSfs(Sfs *sfs){
  int i, n;

  printf("****** SFS ******\n");
  n = sfs->n;
  if(sfs->isFolded)
    n /= 2;
  printf("Polymorphic sites: %d\n", (int)sfs->numPol);
  printf("Monomorphic sites: %d\n", (int)sfs->nullCount);
  printf("n: %d\n", sfs->n);
  printf("u: %.3f\n", sfs->u);
  printf("l: %.3f\n", sfs->l);
    if(sfs->isFolded)
    printf("isFolded: TRUE\n");
  else
    printf("isFolded: FALSE\n");
  if(sfs->type == UNFOLDED)
    printf("type: UNFOLDED\n");
  else if(sfs->type == FOLDED_EVEN)
    printf("type: FOLDED_EVEN\n");
  else if(sfs->type == FOLDED_ODD)
    printf("type: FOLDED_ODD\n");
  else
    printf("type: UNKNOWN\n");
  for(i=0; i<n; i++)
    printf("%d\t%.1f\n", i+1, sfs->f[i]);
  printf("*****************\n");
}

/* splitSfs: Split on site frequency spectrum into args->c site frequency spectra. */
Sfs **splitSfs(Sfs *sfs, gsl_rng *rand, Args *args){
  int *arr;
  int i, j, nc, x, l, c, n;
  Sfs **sfsArr;

  arr = (int *)emalloc((sfs->numPol + sfs->nullCount)*sizeof(int));
  n = sfs->n;
  if(sfs->isFolded)
    n /= 2;
  else
    n--;
  nc = sfs->nullCount;
  for(i=0;i<nc;i++)
    arr[i] = 0;
  x = nc;
  for(i=0;i<n;i++)
    for(j=0;j<sfs->f[i];j++)
      arr[x++] = i+1;
  if(x != sfs->numPol + sfs->nullCount){
    printf("ERROR in splitSfs; sfs->numPol: %d; sfs->nullCount: %d; x: %d; expected(x): %d\n", (int)sfs->numPol, (int)sfs->nullCount, x, (int)(sfs->numPol+sfs->nullCount));
    exit(-1);
  }
  shuffleArr(arr, x, rand);
  l = x / args->c;
  c = 0;
  sfsArr = (Sfs **)emalloc(args->c*sizeof(Sfs *));
  for(i=0;i<args->c;i++){
    sfsArr[i] = newSfs(sfs->n, sfs->type);
    sfsArr[i]->u = sfs->u;
    sfsArr[i]->l = sfs->l;
    for(j=0; j<l; j++){
      if(arr[c]){
	sfsArr[i]->numPol++;
	sfsArr[i]->f[arr[c]-1]++;
      }else
	sfsArr[i]->nullCount++;
      c++;
    }
  }
  free(arr);
  return sfsArr;
}

/* resetSfs: Set the frequencies in sfs to zero. */
void resetSfs(Sfs *sfs){
  int i, n;
  n = sfs->n;
  if(sfs->type == FOLDED_EVEN)
    n /= 2;
  for(i=0; i<n; i++)
    sfs->f[i] = 0;
  sfs->numPol = 0;
  sfs->nullCount = 0;
}

/* addSfs: Add contents of sfs2 to sfs1. */
void addSfs(Sfs *sfs1, Sfs *sfs2){
  int i, n;

  n = sfs1->n;
  if(sfs1->type == FOLDED_EVEN)
    n /= 2;
  sfs1->numPol += sfs2->numPol;
  sfs1->nullCount += sfs2->nullCount;

  for(i=0; i<n; i++)
    sfs1->f[i] += sfs2->f[i];
}

/* normalizeSfs: Normalize SFS such that entries sum to 1. */
void normalizeSfs(Sfs *sfs){
  int i, n;

  n = sfs->n;
  if(sfs->type == FOLDED_EVEN)
    n /= 2;

  for(i=0; i<n; i++)
    sfs->f[i] /= sfs->numPol;
}

/* denormalizeSfs: Reverse normalization. */
void denormalizeSfs(Sfs *sfs){
  int i, n;

  n = sfs->n;
  if(sfs->type == FOLDED_EVEN)
    n /= 2;

  for(i=0; i<n; i++)
    sfs->f[i] *= sfs->numPol;
}
