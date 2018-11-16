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
  Sfs      *s;
  PopSizes *p;
} Rparams;

int newton(Sfs *s, PopSizes *p, Args *a);
double logLik(PopSizes *p, Sfs *s);

#endif
