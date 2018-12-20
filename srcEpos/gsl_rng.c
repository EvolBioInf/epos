/***** gsl_rng.c **********************************
 * Description: Utility functions for random number
 *   generation using the Gnu Scientific Library.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Wed Jun  3 16:39:26 2015
 **************************************************/
#include <stdio.h>
#include <time.h>
#include "gsl_rng.h"

gsl_rng *ini_gsl_rng(Args *args){
  const gsl_rng_type *t;
  gsl_rng *r;
  int idum;
  FILE *fp;

  gsl_rng_env_setup();
  t = gsl_rng_default;
  r = gsl_rng_alloc(t);
  /* seed for random number generation */
  if(args->s != 0){
    idum = args->s;
  }else if((fp = fopen("randomSeed.dat","r")) != NULL){
    if(!fscanf(fp,"%d",&idum))
      printf("WARNING[gsl_rng]: Something went wrong when trying to read the the seed for the random number generator from randomSeed.dat.\n");
    fclose(fp);
  }else
    idum = -time(NULL); 
  gsl_rng_set(r,idum);

  return r;
}


void free_gsl_rng(gsl_rng *r, Args *args){
  FILE *fp;

  if(args->s == 0){
    fp = fopen("randomSeed.dat","w");
    fprintf(fp,"%ld\n",gsl_rng_get(r));
    fclose(fp);
  }
  gsl_rng_free(r);
}
