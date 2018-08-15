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
#include "xval.h"

void scanFile(FILE *fp, Args *args, char *fileName){
  Sfs *sfs, *bSfs;
  PopSizes *ps;
  int i;
  gsl_rng *rand;

  ps = NULL;
  sfs = getSfs(fp, args);
  printSfsStats(sfs);
  rand =  ini_gsl_rng(args);
  ps = getPopSizes(sfs, args);
  printTimes(ps, sfs);
  bSfs = NULL;
  for(i=0; i<args->b; i++){
    printf("Entered bootstrap\n");
    freePopSizes(ps);
    freeSfs(bSfs);
    bSfs = bootstrapSfs(sfs, rand, args);
    printSfsStats(bSfs);
    ps = getPopSizes(bSfs, args);
    printf("#InputFile:\tbootstrapped_%s\n", fileName);
    printTimes(ps, bSfs);
  }
  freeSfs(sfs);
  freeSfs(bSfs);
  freePopSizes(ps);
  free_gsl_rng(rand, args);
}

int main(int argc, char *argv[]){
  int i;
  char *version;
  Args *args;
  FILE *fp;

  version = "0.56";
  setprogname2("epos");
  args = getArgs(argc, argv);
  if(args->v)
    printSplash(version);
  if(args->h || args->e)
    printUsage(version);

  if(args->numInputFiles == 0){
    fp = stdin;
    scanFile(fp, args, "stdin");
  }else{
    for(i=0;i<args->numInputFiles;i++){
      printf("#InputFile:\t%s\n", args->inputFiles[i]);
      fp = efopen(args->inputFiles[i],"r");
      scanFile(fp, args, args->inputFiles[i]);
      fclose(fp);
    }
  }
  free(args);
  free(progname());

  return 0;
}

