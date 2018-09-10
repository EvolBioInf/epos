/***** epos.c *************************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Oct 23 10:10:47 2017
 **************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include "interface.h"
#include "eprintf.h"
#include "sfs.h"
#include "popSizes.h"
#include "util.h"
#include "gsl_rng.h"


void singleAnalysis(Sfs *sfs, Args *args, char *fileName) {
  PopSizes *ps;

  printf("#InputFile:\t");
  if(args->b)
    printf("bootstrapped_");
  printf("%s\n", fileName);
  printSfsStats(sfs);
  ps = getPopSizes(sfs, args);
  printTimes(ps, sfs);
  freePopSizes(ps);
}

void scanFile(FILE *fp, Args *args, char *fileName, gsl_rng *rand){
  Sfs *sfs, *bSfs;

  while((sfs = getSfs(fp, args)) != NULL) {
    if(args->b == 0) {
      singleAnalysis(sfs, args, fileName);
    } else {
      iniBoot(sfs);
      for(int i = 0; i < args->b; i++) {
	bSfs = bootstrapSfs(sfs, rand, args);
	singleAnalysis(bSfs, args, fileName);
	freeSfs(bSfs);
      }
    }
    freeSfs(sfs);
  }
}

int main(int argc, char *argv[]){
  int i;
  char *version;
  Args *args;
  FILE *fp;
  gsl_rng *rand;

  version = "1.0";
  setprogname2("epos");
  args = getArgs(argc, argv);
  if(args->v)
    printSplash(version);
  if(args->h || args->e)
    printUsage(version);

  rand = NULL;
  if(args->b)
    rand = ini_gsl_rng(args);

  if(args->numInputFiles == 0){
    fp = stdin;
    scanFile(fp, args, "stdin", rand);
  }else{
    for(i=0;i<args->numInputFiles;i++){
      fp = efopen(args->inputFiles[i],"r");
      scanFile(fp, args, args->inputFiles[i], rand);
      fclose(fp);
    }
  }
  free(args);
  free(progname());
  if(args->b)
    free_gsl_rng(rand, args);

  return 0;
}

