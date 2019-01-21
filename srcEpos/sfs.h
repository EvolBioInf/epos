/***** sfs.h **************************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Oct 23 10:30:28 2017
 **************************************************/
#ifndef SFS
#define SFS

#include <stdio.h>
#include <gsl/gsl_randist.h>
#include "interface.h"

typedef struct sfs {
  long  *G; /* frequency spectrum          */
  int    a; /* maximum index in G          */
  long   l; /* sequence length             */
  int    n; /* number of haplotypes        */
  double u; /* mutation rate               */
  int    p; /* number of polymorphic sites */
  short  f; /* folded?                     */
} Sfs;

typedef struct sfsSet {
  int n;       /* size of set        */
  Sfs **test;  /* n SFS for testing  */
  Sfs **train; /* n SFS for training */
} SfsSet;

Sfs *readSfs(FILE *fp, Args *args);
void resetReadSfs();
void freeSfs(Sfs *sfs);
void printSfs(Sfs *sfs);
Sfs *newSfs(int n, Args *args);
SfsSet *newSfsSet(Sfs *sfs, Args *args);
void freeSfsSet(SfsSet *ss);
SfsSet *splitSfs(Sfs *sfs, Args *args, gsl_rng *r);

#endif 
