/***** sfs.c **************************************
 * Description: Deal with site frequency spectrum.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Oct 23 10:30:27 2017
 **************************************************/
#include <stdlib.h>
#include <string.h>
#include "eprintf.h"
#include "tab.h"
#include "sfs.h"

short first = 1;
short last  = 0;

void printSfs(Sfs *sfs) {
  printf("****** SFS ******\n");
  printf("Type: ");
  if(sfs->f)
    printf("folded\n");
  else
    printf("unfolded\n");
  printf("Haplotypes: %d\n", sfs->n);
  printf("Sequence length: %ld\n", sfs->l);
  printf("Mutation rate: %g\n", sfs->u);
  printf("# r\tG[r]\n");
  for(int i = 0; i <= sfs->a; i++)
    printf("%d\t%ld\n", i, sfs->G[i]);
  printf("*****************\n");
}

Sfs *newSfs(int n, Args *args) {
  Sfs *sfs = (Sfs *)emalloc(sizeof(Sfs));
  sfs->n = n;
  sfs->G = (long *)emalloc(n * sizeof(long));
  sfs->l = -1;
  sfs->u = args->u;
  sfs->f = 1;
  sfs->p = 0;
  if(args->U)
    sfs->f = 0;
  for(int i = 0; i < n; i++)
    sfs->G[i] = 0.;

  return sfs;
}

void freeSfs(Sfs *sfs){
  if(sfs == NULL)
    return;
  free(sfs->G);
  free(sfs);
}

/* numPol returns the number of polymorphic sites in an SFS */
int numPol(Sfs *sfs) {
  int p = 0;

  for(int i = 1; i <= sfs->a; i++)
    p += sfs->G[i];

  return p;
}

void resetReadSfs() {
  first = 1;
  last  = 0;
}

/* prepSfs prepares the raw data read by readSfs for further analysis */
void prepSfs(Sfs *sfs, int r, Args *args) {
  sfs->u = args->u;
  sfs->l = args->l;
  /* deal with sample size */
  if(sfs->f)
    sfs->n = 2 * r;
  else
    sfs->n = r + 1;
  sfs->a = r;
  sfs->p = numPol(sfs);
  /* deal with monomorphic sites */
  if(sfs->G[0] == 0) {
    if(sfs->l == 0) {
      fprintf(stderr, "ERROR[epos]: Please include either the zero-clas in the SFS or enter the sequence length via the -l option.\n");
      exit(-1);
    } 
    sfs->G[0] = args->l - numPol(sfs);
    if(sfs->G[0] <= 0) {
      fprintf(stderr, "ERROR[epos]: The sequence length (%ld) must be greater than the number of polymorphic sites (%d).\n", sfs->l, numPol(sfs));
      exit(-1);
    }
  } else
    sfs->l = sfs->G[0] + numPol(sfs);
}

/* getSfs obtains the next SFS from an open file */
Sfs *readSfs(FILE *fp, Args *args){
  char *line;
  int r = 0;
  double f;

  if(last)
    return NULL;
  Sfs *sfs = newSfs(1, args);
  while((line = tabGetLine(fp)) != NULL){
    if(line[0] == '#') {          /* header line */
      if(first) {
	first = 0;
	continue;
      } else {
	prepSfs(sfs, r, args);
	return sfs;
      }
    }
    r = atoi(tabField(0));       /* degree */
    f = atof(tabField(1));       /* frequency */
    sfs->G = (long *)erealloc(sfs->G, (r + 1) * sizeof(long));
    sfs->G[r] = f;
  }
  if(line == NULL)
    last = 1;
  prepSfs(sfs, r, args);

  return sfs;
}
