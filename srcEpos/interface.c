/***** interface.c ************************************************
 * Description: Routine for gathering arguments from the command
 *              line.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Sun Jun 20 13:12:10 2004.
 *****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "interface.h"
#include "eprintf.h"

Args *newArgs() {
  Args *args = (Args *)emalloc(sizeof(Args));
  args->h = 0;
  args->v = 0;
  args->e = 0;
  args->U = 0;
  args->t = 0;
  args->u = DEFAULT_U;
  args->c = -1;
  args->l = 0;
  args->E = DEFAULT_E;
  args->x = DEFAULT_X;
  args->s = 0;
  args->L = NULL;
  args->al = (int *)emalloc(3 * sizeof(int));
  args->nl = 0;
  return args;
}

void freeArgs(Args *args) {
  if(args->L != NULL)
    free(args->L);
  free(args->al);
  free(args);
}

void extractLevels(Args *args) {
  char *c;
  c = strtok(args->L, ",");

  int i = atoi(c);
  if(i != 2) {
    args->nl++;
    args->al[args->nl] = 2;
  }
  args->nl++;
  args->al[args->nl] = i;
  while((c = strtok(NULL, ",")) != NULL){
    i = atoi(c);
    args->nl++;
    args->al = (int *)erealloc(args->al, (args->nl + 1) * sizeof(int));
    args->al[args->nl] = i;
  }
}

Args *getArgs(int argc, char *argv[]){
  int c;
  char *optString = "hvUtu:l:L:c:E:x:s:";

  Args *args = newArgs();
  c = getopt(argc, argv, optString);
  while(c != -1){
    switch(c){
      break;
    case 'l':                           /* sequence length */
      args->l = (long)atof(optarg);
      break;
    case 'L':                           /* preset levels */
      args->L = estrdup(optarg);
      break;
    case 'E':                           /* levels of exhaustive search */
      args->E = atoi(optarg);
      break;
    case 's':                           /* seed for random number generator */
      args->s = atoi(optarg);
      break;
    case 'x':                           /* number of categories for cross validation */
      args->x = atoi(optarg);
      if(args->x < 1)
	args->x = 1;
      break;
    case 'U':                           /* unfolded */
      args->U = 1;
      break;
    case 't':                           /* testing? */
      args->t = 1;
      break;
    case 'u':
      args->u = atof(optarg);           /* mutation rate */
      break;
    case 'c':
      args->c = atof(optarg);           /* minimum logLik change for acceptance of new level */
      if(args->c < 0)
	args->c = 0.;
      break;
    case '?':                           /* fall-through is intentional */
    case 'h':                           /* print help */
      args->h = 1;
      break;
    case 'v':                           /* print version */
      args->v = 1;
      break;
    default:
      printf("# unknown argument: %c\n",c);
      args->e = 1;
      return args;
    }
    c = getopt(argc, argv, optString);
  }
  args->inputFiles = argv + optind;
  args->numInputFiles = argc - optind;
  if(args->L)
    extractLevels(args);
  if(args->c < 0) {
    if(args->x == 1)
      args->c = DEFAULT_C;
    else
      args->c = 0.;
  }
  return args;
}


void printUsage(){
  printf("Usage: %s [options] [inputFiles]\n",progname());
  printf("Estimate population size from site frequency data\n");
  printf("Example: ./epos -l 10000000 -u 1.2e-8 data/testNewtonF.dat\n");
  printf("Options:\n");
  printf("\t[-l NUM sequence length; default: include zero-class in SFS]\n");
  printf("\t[-u NUM per nucleotide mutation rate; default: %g]\n", DEFAULT_U);
  printf("\t[-c NUM minimum change in log-likelihood for acceptance of new level; default: %g if -x 1, 0 otherwise]\n", DEFAULT_C);
  printf("\t[-E NUM levels searched exhaustively; default: greedy search; -E > 2 differs from greedy]\n");
  printf("\t[-L NUM1,NUM2,... use preset levels NUM1,NUM2,...; default: search for optimal levels]\n");
  printf("\t[-x NUM categories for cross validation; default: %d]\n", DEFAULT_X);
  printf("\t[-s NUM seed for random number generator; default: system]\n");
  printf("\t[-U unfolded site frequency spectrum; default: folded]\n");
  printf("\t[-t exectute test routines for debuggin]\n");
  printf("\t[-h print this help message and exit]\n");
  printf("\t[-v print program information and exit]\n");
  exit(0);
}

void printSplash(char *version, char *date){
  printf("%s ", progname());
  int l = strlen(version);
  for(int i = 0; i < l - 1; i++)
    printf("%c", version[i]);
  printf(", %s\n", date);
  printf("Written by Bernhard Haubold.\n");
  printf("Distributed under the GNU General Public License.\n");
  printf("Please send bug reports to haubold@evolbio.mpg.de\n");
  exit(0);
}

