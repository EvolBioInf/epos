/***** util.h *************************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Fri Dec 15 10:36:58 2017
 **************************************************/
#ifndef UTIL
#define UTIL

#include "popSizes.h"
#include "sfs.h"
#include "interface.h"

inline double binomial(int n, int k);
double fi(int i, int n, int r, int *k);
double gi(int i, int n, int r, int *k);
double hi(int i, int n, int r, int *k);
int delta(int j, int k);
void printTimes(PopSizes *ps, Sfs *sfs);
void printAaa(PopSizes *ps, Sfs *sfs);
void printSfsStats(Sfs *sfs);
double watterson(Sfs *sfs, Args *args);
void iniBinom(int n);
void freeBinom();
inline int max(int a, int b);
void shuffle(int *a, long n, gsl_rng *r);
void testShuffle(int n, gsl_rng *r);

extern double **bin;

inline double binomial(int n, int k){

  if(__builtin_expect(n < k, 0))
    return 0;
  else
    return bin[n][k];
}

inline int max(int a, int b) {
  if(a > b)
    return a;
  else
    return b;
}

#endif
