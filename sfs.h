/***** sfs.h **************************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Oct 23 10:30:28 2017
 **************************************************/
#ifndef SFS
#define SFS
#define UNFOLDED 0
#define FOLDED_EVEN 1
#define FOLDED_ODD 2

#include <stdio.h>
#include "interface.h"
#include "gsl_rng.h"

/* define site frequency spectrum */
typedef struct sfs{
  double *f;        /* frequency spectrum */
  double numPol;    /* number of polymorphic sites */
  int n;            /* sample size */
  double nullCount; /* number of unmutated sites */
  double u;         /* mutation rate */
  double l;         /* lambda */
  double iniP;      /* initial population size */
  short type;       /* type of site frequency spectrum */
  char isFolded;
} Sfs;

Sfs *getSfs(FILE *fp, Args *args);
void freeSfs(Sfs *sfs);
Sfs *bootstrapSfs(Sfs *sfs, gsl_rng *rand, Args *args);
void printSfs(Sfs *sfs);
Sfs **splitSfs(Sfs *sfs, gsl_rng *rand, Args *args);
void resetSfs(Sfs *sfs);
void addSfs(Sfs *sfs1, Sfs *sfs2);
Sfs *newSfs(int n, int type);
void normalizeSfs(Sfs *sfs);
void denormalizeSfs(Sfs *sfs);

#endif 
