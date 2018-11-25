/***** epos.c *************************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Oct 23 10:10:47 2017
 **************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include "tab.h"
#include "interface.h"
#include "eprintf.h"
#include "sfs.h"
#include "popSizes.h"
#include "util.h"
#include "gsl_rng.h"
#include "search.h"

void test(Args *args);

void analysis(Sfs *sfs, Args *args, char *fileName) {
  PopSizes *ps;

  printf("#InputFile:\t");
  printf("%s\n", fileName);
  /* printSfs(sfs); */
  printSfsStats(sfs);
  ps = searchLevels(sfs, args);
  printTimes(ps, sfs);
  /* if(args->a) */
  /*   printAaa(ps, sfs); */
  /* else */
  /*   printTimes(ps, sfs); */
  /* freePopSizes(ps); */
}

void scanFile(FILE *fp, Args *args, char *fileName){
  Sfs *sfs;

  while((sfs = readSfs(fp, args)) != NULL) {
    iniBinom(sfs->n);
    analysis(sfs, args, fileName);
    freeBinom();
    freeSfs(sfs);
  }
  tabFree();
}

int main(int argc, char *argv[]){
  char *version;
  Args *args;
  FILE *fp;

  args = getArgs(argc, argv);
  if(args->t) {
    test(args);
    return 0;
  }
  version = "1.0";
  setprogname2("epos");
  if(args->v)
    printSplash(version);
  if(args->h || args->e)
    printUsage(version);
  if(args->numInputFiles == 0){
    fp = stdin;
    scanFile(fp, args, "stdin");
  }else{
    for(int i = 0; i < args->numInputFiles; i++){
      fp = efopen(args->inputFiles[i], "r");
      scanFile(fp, args, args->inputFiles[i]);
      fclose(fp);
    }
  }
  free(args);
  free(progname());

  return 0;
}

