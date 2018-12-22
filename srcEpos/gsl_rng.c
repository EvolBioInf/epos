/***** gsl_rng.c **********************************
 * Description: Utility functions for random number
 *   generation using the Gnu Scientific Library.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Wed Jun  3 16:39:26 2015
 **************************************************/
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include "gsl_rng.h"

gsl_rng *ini_gsl_rng(Args *args){
  const gsl_rng_type *t;
  gsl_rng *r;
  int seed;

  gsl_rng_env_setup();
  t = gsl_rng_default;
  r = gsl_rng_alloc(t);
  /* seed for random number generation */
  if(args->s != 0)
    seed = args->s;
  else
    seed = (unsigned int)((unsigned) time(NULL) + getpid());
  gsl_rng_set(r, seed);

  return r;
}
