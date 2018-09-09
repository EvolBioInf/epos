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

void scanFile(FILE *fp, Args *args, char *fileName){
  Sfs *sfs;
  PopSizes *ps;

  ps = NULL;
  while((sfs = getSfs(fp, args)) != NULL) {
    printSfsStats(sfs);
    ps = getPopSizes(sfs, args);
    printTimes(ps, sfs);
    freeSfs(sfs);
    freePopSizes(ps);
  }
}

int main(int argc, char *argv[]){
  int i;
  char *version;
  Args *args;
  FILE *fp;

  version = "1.0";
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

