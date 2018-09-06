/***** newton.h ***********************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Tue Jul 24 18:12:30 2018
 **************************************************/
#ifndef NEWTON
#define NEWTON

#include "interface.h"
#include "sfs.h"
#include "popSizes.h"

typedef struct rparams {
  int n;       /* sample size                          */
  int m;       /* number of population sizes           */
  int l;       /* sequence length                      */
  int o;       /* number of unmutated positions        */
  int *k;      /* locations of population size changes */
  double *g;   /* observed site frequency spectrum     */
  double u;    /* mutation rate                        */
} Rparams;

int newton(Sfs *sfs, PopSizes *ps);
double logLik(PopSizes *ps, Sfs *sfs);

#endif
