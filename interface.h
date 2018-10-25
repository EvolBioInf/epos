#ifndef INTERFACE
#define INTERFACE

#include <float.h>

#define DEFAULT_U 5.e-9 /* default mutation rate */
#define DEFAULT_C 2.    /* default minimum change in log-likelihood for acceptance of new level */ 
/* define argument container */
typedef struct args{
  char h;      /* help message? */
  char v;      /* version message? */
  char e;      /* error message? */
  char U;      /* unfolded? */
  char a;      /* average age of an allele? */
  int l;       /* sequence length */
  int s;       /* seed for random number generator */
  int b;       /* number of bootstrap replicates */
  char **inputFiles;
  char *x;     /* excluded frequency categories, e. g. 1,2 for excluding singletons & doubletons */
  double u;    /* mutation rate */
  double c;    /* minimum change */
  int numInputFiles;
} Args;

Args *getArgs(int argc, char *argv[]);
void printUsage(char *version);
void printSplash(char *version);
void freeArgs();

#endif
