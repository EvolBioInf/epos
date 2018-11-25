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
  int     m;  /* number of population size changes, 1 <= m <= n-1 */
  double *N;  /* m population sizes                               */
  int    *k;  /* m+1 positions where population size changes      */
  double  l;  /* log-likelihood                                   */
} PopSizes;

PopSizes *newPopSizes(Sfs *sfs);
void freePopSizes(PopSizes *ps);
double logLik(PopSizes *p, Sfs *s);
double expF(PopSizes *p, Sfs *s, int r);
double expG(PopSizes *p, Sfs *s, int r);
int negPopSizes(PopSizes *ps);

#endif
