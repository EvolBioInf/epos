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
  args->h  = 0;
  args->v  = 0;
  args->e  = 0;
  args->U  = 0;
  args->t  = 0;
  args->u  = DEFAULT_U;
  args->c  = -1;
  args->l  = 0;
  args->E  = DEFAULT_E;
  args->m  = 0;
  args->k  = DEFAULT_K;
  args->s  = 0;
  args->o  = 0;
  args->d  = 0;
  args->L  = NULL;
  args->al = NULL;
  args->nl = 0;
  args->x  = NULL;
  args->ax = NULL;
  args->nx = 0;

  return args;
}

void freeArgs(Args *args) {
  if(args->L != NULL)
    free(args->L);
  if(args->al != NULL)
    free(args->al);
  if(args->x != NULL)
    free(args->x);
  if(args->ax != NULL)
    free(args->ax);
  free(args);
}

void freeInts(Ints *in) {
  free(in->a);
  free(in);
}

Ints *newInts() {
  Ints *in = emalloc(sizeof(Ints));
  in->a = NULL;
  in->n = 0;

  return in;
}

/* cmpInt compares integers */
int cmpInt (const void *a, const void *b) {
  const int *i1 = (const int *)a;
  const int *i2 = (const int *)b;
  return *i1 - *i2;
}


Ints *extractInts(char *s) {
  Ints *in =  newInts();
  char *c = strtok(s, ",");
  in->a = (int *)emalloc(sizeof(int));
  in->a[in->n++] = atoi(c);
  while((c = strtok(NULL, ",")) != NULL){
    in->a = (int *)erealloc(in->a, (in->n + 1) * sizeof(int));
    in->a[in->n++] = atoi(c);
  }
  qsort(in->a, in->n, sizeof(int), cmpInt);
  return in;
}

void extractClasses(Args *args) {
  Ints *in = extractInts(args->x);
  args->ax = (int *)malloc((in->n + 1) * sizeof(int));
  args->nx = in->n;

  for(int i = 0; i < in->n; i++)
    args->ax[i] = in->a[i];

  freeInts(in);
}

void extractLevels(Args *args) {
  Ints *in = extractInts(args->L);
  args->al = (int *)malloc((in->n + 1) * sizeof(int));
  args->nl = in->n;
  
  for(int i = 0; i < in->n; i++) {
    args->al[i + 1] = in->a[i];
    printf("al[%d]: %d\n", i + 1, args->al[i + 1]);
  }

  freeInts(in);
}


Args *getArgs(int argc, char *argv[]){
  int c;
  char *optString = "hvUtdou:l:L:c:E:x:s:k:m:";

  Args *args = newArgs();
  c = getopt(argc, argv, optString);
  while(c != -1){
    switch(c){
      break;
    case 'l':                           /* sequence length */
      args->l = (long)atof(optarg);
      break;
    case 'm':                           /* maximal level searched exhaustively */
      args->m = atoi(optarg);
      break;
    case 'L':                           /* preset levels */
      args->L = estrdup(optarg);
      extractLevels(args);
      break;
    case 'o':                           /* print obs./exp. freq. spec? */
      args->o = 1;
      break;
    case 'd':                           /* print debug information */
      args->d = 1;
      break;
    case 'E':                           /* levels of exhaustive search */
      args->E = atoi(optarg);
      break;
    case 's':                           /* seed for random number generator */
      args->s = atoi(optarg);
      break;
    case 'k':                           /* number of categories for cross validation */
      args->k = atoi(optarg);
      if(args->k < 1)
	args->k = 1;
      break;
    case 'x':                           /* excluded frequencies */
      args->x = estrdup(optarg);
      extractClasses(args);
      for(int i = 0; i < args->nx; i++)
	if(args->ax[i] < 0)
	  printf("getArgs - 1\n");
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

  if(args->c < 0) {
    if(args->k == 1)
      args->c = DEFAULT_C;
    else
      args->c = 0.;
  }
  return args;
}


void printUsage(){
  printf("Usage: %s [options] [inputFiles]\n",progname());
  printf("Estimate population size from site frequency data\n");
  printf("Example: ./epos -l 1e7 -u 1.2e-8 data/testNewtonF.dat\n");
  printf("Options:\n");
  printf("\t[-l NUM sequence length; default: include zero-class in SFS]\n");
  printf("\t[-u NUM per nucleotide mutation rate; default: %g]\n", DEFAULT_U);
  printf("\t[-c NUM minimum change in log-likelihood for acceptance of new level; default: %g if -k 1, 0 otherwise]\n", DEFAULT_C);
  printf("\t[-E NUM levels searched exhaustively; default: greedy search; -E > 2 differs from greedy]\n");
  printf("\t[-m NUM maximal level searched exhaustively; default sample size]\n");
  printf("\t[-L NUM1,NUM2,... use preset levels NUM1,NUM2,...; default: search for optimal levels]\n");
  printf("\t[-x NUM1,NUM2,... exclude NUM1-ers, NUM2-ers; default: include all frequency classes]\n");
  printf("\t[-k NUM cross validation with NUM categories: no cross-validation]\n");
  printf("\t[-s NUM seed for random number generator; default: system]\n");
  printf("\t[-U unfolded site frequency spectrum; default: folded]\n");
  printf("\t[-o print observed and expected site frequency spectrum]\n");
  printf("\t[-d print debug information: split site frequency spectra]\n");
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

