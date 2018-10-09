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

Args *args;

Args *getArgs(int argc, char *argv[]){
  char c;
  char *optString = "hvUau:l:c:s:b:x:";

  args = (Args *)emalloc(sizeof(Args));
  args->h = 0;
  args->v = 0;
  args->e = 0;
  args->U = 0;
  args->s = 0;
  args->b = 0;
  args->a = 0;
  args->x = NULL;
  args->u = DEFAULT_U;
  args->c = DEFAULT_C;
  args->l = 0;

  c = getopt(argc, argv, optString);
  while(c != -1){
    switch(c){
      break;
    case 'l':                           /* sequence length */
      args->l = atoi(optarg);
      break;
    case 's':                           /* seed for random number generator */
      args->s = atoi(optarg);
      break;
    case 'b':                           /* number of bootstrap replicates */
      args->b = atoi(optarg);
      break;
    case 'U':                           /* unfolded */
      args->U = 1;
      break;
    case 'x':
      args->x = estrdup(optarg);
      break;
    case 'a':                           /* average age of an allele */
      args->a = 1;
      break;
    case 'u':
      args->u = atof(optarg);           /* mutation rate */
      break;
    case 'c':
      args->c = atof(optarg);           /* minimum logLik change for acceptance of new level */
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
  return args;
}


void printUsage(char *version){
  printf("Usage: %s [options] [inputFiles]\n",progname());
  printf("Estimate population size from site frequency data\n");
  printf("Usage: ./epos [options] sfs.txt\n");
  printf("Example: ./epos -l 10000000 -u 1.2e-8 data/testNewtonF.dat\n");
  printf("Options:\n");
  printf("\t[-l NUM sequence length; default: include zero-class in SFS]\n");
  printf("\t[-u NUM per nucleotide mutation rate; default: %g]\n", DEFAULT_U);
  printf("\t[-c NUM minimum change in log-likelihood for accepatnce of new level; default: %g]\n", DEFAULT_C);
  printf("\t[-b NUM number of bootstrap replicates; default: no bootrstrap]\n");
  printf("\t[-s NUM seed for random number generator used in bootstrap; default: system, randomSeed.dat]\n");
  printf("\t[-x NUM1,NUM2... exclude frequency categories NUM1, NUM2, etc.; default: include all categories]\n");
  printf("\t[-U unfolded site frequency spectrum; default: folded]\n");
  printf("\t[-a print average ages of alleles; default: print population size]\n");
  printf("\t[-h print this help message and exit]\n");
  printf("\t[-v print program information and exit]\n");
  exit(0);
}

void printSplash(char *version){
  printf("%s %s\n",progname(),version);
  printf("Written by Bernhard Haubold.\n");
  printf("Distributed under the GNU General Public License.\n");
  printf("Please send bug reports to haubold@evolbio.mpg.de\n");
  exit(0);
}

