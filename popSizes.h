/***** popSizes.h *********************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Dec 18 11:19:35 2017
 **************************************************/
#ifndef POPSIZES
#define POPSIZES

#include "sfs.h"
#include "interface.h"

typedef struct popSizes {
  double *N;
  double *allN;
  double psi;
  double watterson;
  int m;
  int *k;
  int *prevK;
  double *iniN; /* initial population sizes */
  double *aaa;  /* average ages of alleles */
  double *asa;  /* average sizes of alleles */
  int prevM;
  int n;
} PopSizes;

PopSizes *getPopSizes(Sfs *sfs, Args *args);
PopSizes *newPopSizes(Sfs *sfs);
void addTestK(PopSizes *ps, int k);
void addK(PopSizes *ps, int k);
void restoreK(PopSizes *ps);
void freePopSizes(PopSizes *ps);
int negPopSizes(PopSizes *ps);
double psi(PopSizes *ps, Sfs *sfs);
void findIniN(PopSizes *ps, Sfs *sfs);
void compAaa(PopSizes *ps, Sfs *sfs);

#endif
