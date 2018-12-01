/***** interface.c ************************************************
 * Description: Routine for gathering arguments from the command
 *              line.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Sun Jun 20 13:12:10 2004.
 *****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "interface.h"
#include "eprintf.h"

Args *args;

Args *getArgs(int argc, char *argv[]){
  int c;
  char *optString = "hvn:";

  args = (Args *)emalloc(sizeof(Args));
  args->n = 0;
  args->h = 0;
  args->v = 0;
  args->e = 0;

  c = getopt(argc, argv, optString);
  while(c != -1){
    switch(c){
    case 'n':                           /* sample size */
      args->n = atoi(optarg);
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
  if(!args->n)
    args->e = 1;
  return args;
}


void printUsage(char *version){
  printf("Usage: %s [options] [inputFiles]\n",progname());
  printf("Compute the average ages of alleles and the average population size they experienced.\n");
  printf("Example: epos kap144i.dat | epos2ages -n 144\n");
  printf("Options:\n");
  printf("\t-n <NUM> sample size (number of haplotypes surveyed)\n");
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
