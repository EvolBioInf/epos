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
  char *optString = "nhvVUNc:u:s:l:";

  args = (Args *)emalloc(sizeof(Args));
  args->h = 0;
  args->v = 0;
  args->V = 0;
  args->e = 0;
  args->b = 0;
  args->U = 0;
  args->s = 0;
  args->m = 0;
  args->n = 0; /* always allow negative population sizes */
  args->N = 0;
  args->c = DEFAULT_C;
  args->L = 0;
  args->f = DEFAULT_F;
  args->d = DEFAULT_D;
  args->u = DEFAULT_U;
  args->p = 0;
  args->l = DEFAULT_L;

  c = getopt(argc, argv, optString);
  while(c != -1){
    switch(c){
      break;
    case 's':
      args->s = atoi(optarg);           /* seed for random number generator */
      break;
    case 'f':                           /* factor for scaling lambda */
      args->f = atof(optarg);
      break;
    case 'c':                           /* number of categories for cross-validation */
      args->c = atoi(optarg);
      break;
    case 'p':                           /* print matrix */
      args->p = 1;
      break;
    case 'n':                           /* allow negative population sizes */
      args->n = 1;
      break;
    case 'N':                           /* Newton procedure */
      args->N = 1;
      break;
    case 'L':                           /* search for optimal lambda */
      args->L = 1;
      break;
    case 'V':                           /* verbose */
      args->V = 1;
      break;
    case 'U':                           /* unfolded */
      args->U = 1;
      break;
    case 'u':
      args->u = atof(optarg);           /* mutation rate */
      break;
    case 'l':
      args->l = atof(optarg);           /* lambda */
      break;
    case 'b':
      args->b = atoi(optarg);           /* number of bootstrap replicates */
      break;
    case 'm':                           /* maximum number of population sizes */
      args->m = atoi(optarg); 
      break;
    case 'd':
      args->d = atof(optarg);           /* minimal decrease in objective function */
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
  if(args->c == 1)   /* No cross-validation */
    args->d = CHI_THRESHOLD;
  return args;
}


void printUsage(char *version){
  printf("Usage: %s [options] [inputFiles]\n",progname());
  printf("Estimate population size from site frequency data\n");
  printf("Example: epos sfs.txt\n");
  printf("Options:\n");
  printf("\t[-u NUM per nucleotide mutation rate; default: %g; estimated from SFS if zero-class present]\n", DEFAULT_U);
  printf("\t[-l NUM lambda; default: %.3g]\n",DEFAULT_L);
  printf("\t[-c NUM number of categories for cross-validation; default: %d]\n", DEFAULT_C);
  printf("\t[-s NUM seed for random number generator; default: time, file]\n");
  printf("\t[-n allow negative population sizes]\n");
  printf("\t[-N Newton procedure; default: linear algebra]\n");
  printf("\t[-U unfolded site frequency spectrum as input]\n");
  printf("\t[-V verbose output for debugging]\n");
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

