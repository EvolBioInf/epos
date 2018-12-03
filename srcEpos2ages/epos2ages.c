/***** epos2ages.c ********************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Thu Nov 29 12:21:20 2018
 **************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "interface.h"
#include "eprintf.h"
#include "popSizes.h"
#include "ages.h"
#include "util.h"
#include "tab.h"

void scanFile(FILE *fp, Args *args){
  PopSizes *p;
  Ages *a;

  while((p = nextPs(fp)) != NULL) {
    p->k[p->m + 1] = args->n + 1;
    a = compAges(p);
    printAges(a);
    freeAges(a);
    freePopSizes(p);
  }
}

int main(int argc, char *argv[]){
  int i;
  Args *args;
  FILE *fp;

  char *version = VERSION;
  char *date    = DATE;
  setprogname2("epos2age");
  args = getArgs(argc, argv);
  iniBinom(args->n);
  if(args->v)
    printSplash(version, date);
  if(args->h || args->e)
    printUsage(version);
  if(args->numInputFiles == 0){
    fp = stdin;
    scanFile(fp, args);
  }else{
    for(i=0;i<args->numInputFiles;i++){
      fp = efopen(args->inputFiles[i],"r");
      scanFile(fp, args);
      fclose(fp);
    }
  }
  free(args);
  free(progname());
  freeBinom();
  tabFree();
  return 0;
}

