#ifndef INTERFACE
#define INTERFACE

#include <float.h>

#define DEFAULT_U 5.e-9         /* default mutation rate */
/* define argument container */
typedef struct args{
  char h;      /* help message? */
  char v;      /* version message? */
  char e;      /* error message? */
  char U;      /* unfolded? */
  char w;      /* watterson's estimator is initial population size */
  int l;       /* sequence length */
  char **inputFiles;
  double u;    /* mutation rate */
  int numInputFiles;
} Args;

Args *getArgs(int argc, char *argv[]);
void printUsage(char *version);
void printSplash(char *version);

#endif
