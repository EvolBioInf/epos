/***** testNewton.c *******************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Fri Nov 23 11:30:37 2018
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
#include "newton.h"

void scanFile(FILE *fp, Args *args, char *fileName){
  Sfs *sfs;
  PopSizes *ps;

  sfs = readSfs(fp, args);
  iniBinom(sfs->n);
  ps = newPopSizes(sfs);
  ps->m = 4;
  ps->k[1] = 2;
  ps->k[2] = 4;
  ps->k[3] = 6;
  ps->k[4] = 25;
  ps->k[ps->m + 1] = sfs->n + 1;
  newton(sfs, ps, args);
  printSfsStats(sfs);
  printTimes(ps, sfs);
  freeBinom();
  freeSfs(sfs);
  tabFree();
}

int main(int argc, char *argv[]){
  char *version;
  Args *args;
  FILE *fp;

  version = "0.1";
  setprogname2("testEpos");
  args = getArgs(argc, argv);
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

