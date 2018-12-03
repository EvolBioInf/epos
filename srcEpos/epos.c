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
#include "search.h"
#include "newton.h"

void test(Args *args);

PopSizes *presetLevels(Sfs *sfs, Args *args) {
    PopSizes *ps = newPopSizes(sfs);
    for(int i = 1; i <= args->nl; i++)
      ps->k[i] = args->al[i];
    ps->k[args->nl + 1] = sfs->n + 1;
    ps->m = args->nl;
    newton(sfs, ps, args);
    return ps;
}

void analysis(Sfs *sfs, Args *args, char *fileName) {
  PopSizes *ps;

  printf("#InputFile:\t");
  printf("%s\n", fileName);
  printSfsStats(sfs);
  if(args->L)
    ps = presetLevels(sfs, args);
  else {
    ps = searchLevels(sfs, args);
  }
  printTimes(ps, sfs);
  freePopSizes(ps);
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
  Args *args;
  FILE *fp;

  args = getArgs(argc, argv);
  if(args->t) {
    test(args);
    return 0;
  }
  char *version = VERSION;
  char *date    = DATE;
  setprogname2("epos");
  if(args->v)
    printSplash(version, date);
  if(args->h || args->e)
    printUsage(version, date);
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
  freeArgs(args);
  free(progname());

  return 0;
}

