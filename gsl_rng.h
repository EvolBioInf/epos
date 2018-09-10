/***** gsl_rng.h **********************************
 * Description: Utility functions for random number
 *   generation using the Gnu Scientific Library.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Wed Jun  3 16:39:24 2015
 **************************************************/
#ifndef GSL_RNG
#define GSL_RNG

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include "interface.h"

gsl_rng *ini_gsl_rng(Args *args);
void free_gsl_rng(gsl_rng *r, Args *args);

#endif
