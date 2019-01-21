/***** util.h *************************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Fri Dec 15 10:36:58 2017
 **************************************************/
#ifndef UTIL
#define UTIL

#include "popSizes.h"
#include "sfs.h"

double binomial(int n, int k);
double fi(int i, int n, int r, int *k);
double gi(int i, int n, int r, int *k);
double hi(int i, int n, int r, int *k);
int delta(int j, int k);
void printTimes(PopSizes *ps, Sfs *sfs);
void printAaa(PopSizes *ps, Sfs *sfs);
void printSfsStats(Sfs *sfs);
double watterson(Sfs *sfs);
void iniBinom(int n);
void freeBinom();
int max(int a, int b);
void shuffle(int *a, long n, gsl_rng *r);
void testShuffle(int n, gsl_rng *r);

#endif
