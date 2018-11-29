/***** popSizes.h *********************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Thu Nov 29 12:22:28 2018
 **************************************************/
#ifndef POPSIZES
#define POPSIZES

#include <stdio.h>

typedef struct popSizes {
  int     m;  /* number of population size changes, 1 <= m <= n-1 */
  double *N;  /* m population sizes                               */
  int    *k;  /* m+1 positions where population size changes      */
} PopSizes;

PopSizes *nextPs(FILE *fp);
void freePopSizes(PopSizes *p);

#endif
