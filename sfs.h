/***** sfs.h **************************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Oct 23 10:30:28 2017
 **************************************************/
#ifndef SFS
#define SFS

#include <stdio.h>
#include "interface.h"

typedef struct sfs{
  int   *G; /* frequency spectrum          */
  int    a; /* maximum index in G          */
  int    l; /* sequence length             */
  int    n; /* number of haplotypes        */
  double u; /* mutation rate               */
  int    p; /* number of polymorphic sites */
  short  f; /* folded?                     */
} Sfs;

Sfs *readSfs(FILE *fp, Args *args);
void freeSfs(Sfs *sfs);
void printSfs(Sfs *sfs);
Sfs *newSfs(int n, Args *args);

#endif 
