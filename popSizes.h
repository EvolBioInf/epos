/***** popSizes.h *********************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Dec 18 11:19:35 2017
 **************************************************/
#ifndef POPSIZES
#define POPSIZES

#include "sfs.h"
#include "interface.h"

typedef struct popSizes{
  double *N;
  double psi;
  int m;
  int *k;
  int *prevK;
  int prevM;
  int n;
}PopSizes;

PopSizes *getPopSizes(Sfs *sfs, Args *args);
PopSizes *newPopSizes(Sfs *sfs);
void addTestK(PopSizes *ps, int k);
void addK(PopSizes *ps, int k);
void restoreK(PopSizes *ps);
void freePopSizes(PopSizes *ps);
int negPopSizes(PopSizes *ps);
double psi(PopSizes *ps, Sfs *sfs);

#endif
