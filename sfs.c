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
    sfs->f[i] = 0;
  sfs->type = type;
  sfs->nullCount = 0.;
  sfs->numPol = 0.;
  sfs->u = 0.;
  sfs->iniP = 0;

  return sfs;
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

  if(args->U){
    type = UNFOLDED;
  }else{
    type = FOLDED_EVEN;
  }
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
