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

/* define site frequency spectrum */
typedef struct sfs{
  double *f;        /* frequency spectrum */
  double numPol;    /* number of polymorphic sites */
  int n;            /* sample size */
  double nullCount; /* number of unmutated sites */
  double u;         /* mutation rate */
  short type;       /* type of site frequency spectrum */
} Sfs;

Sfs *getSfs(FILE *fp, Args *args);
void freeSfs(Sfs *sfs);
void printSfs(Sfs *sfs);
void resetSfs(Sfs *sfs);
Sfs *newSfs(int n, int type);

#endif 
