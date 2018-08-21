#ifndef INTERFACE
#define INTERFACE

#include <float.h>

#define DEFAULT_U 5.e-9         /* default mutation rate */
#define DEFAULT_D DBL_EPSILON   /* default delta */
#define DEFAULT_F 1.0           /* default factor for lambda scaling */
#define DEFAULT_L 0.0           /* default lambda */
#define DEFAULT_C 5             /* default number of categories for cross-validation */
#define CHI_THRESHOLD 3.84      /* P=0.05, 1 degree of freedom */
/* define argument container */
typedef struct args{
  char h;      /* help message? */
  char v;      /* version message? */
  char e;      /* error message? */
  int E;       /* sequence length */
  char p;      /* print matrix? */
  char V;      /* verbose output? */
  char U;      /* unfolded? */
  char n;      /* allow negative population sizes */
  char N;      /* Newton procedure? */
  char L;      /* search lambda through cross-validation */
  int b;       /* bootstrap */
  int s;       /* seed for random number generator */
  int m;       /* maximum number of population sizes */
  int c;       /* number of categories for cross-validation */
  char **inputFiles;
  double u;    /* mutation rate */
  double l;    /* lambda */
  double d;    /* delta for estimation */
  double f;    /* factor for lambda scaling */
  int numInputFiles;
} Args;

Args *getArgs(int argc, char *argv[]);
void printUsage(char *version);
void printSplash(char *version);

#endif
