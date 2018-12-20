/***** search.h ***********************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Sat Nov 24 14:39:42 2018
 **************************************************/
#ifndef SEARCH
#define SEARCH

#include <gsl/gsl_rng.h>

PopSizes *searchLevels(Sfs *sfs, Args *args, gsl_rng *r);

#endif
