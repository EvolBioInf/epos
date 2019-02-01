#ifndef INTERFACE
#define INTERFACE

#include <float.h>

#define DEFAULT_U 5.e-9 /* mutation rate                                                */
#define DEFAULT_C 2.    /* minimum change in log-likelihood for acceptance of new level */
#define DEFAULT_E 2     /* number of levels searched exhaustively                       */
#define DEFAULT_X 5     /* number of categories for cross-validation                    */
/* define argument container */
typedef struct args{
  char h;      /* help message?                             */
  char v;      /* version message?                          */
  char e;      /* error message?                            */
  char U;      /* unfolded?                                 */
  char t;      /* execute test routines?                    */
  char o;      /* print observed/expected freq. spec?       */
  char d;      /* print debug information?                  */
  long l;      /* sequence length                           */
  char *L;     /* preset levels                             */
  int nl;      /* number of levels                          */
  int *al;     /* array of nl levels                        */
  char *X;     /* excluded frequency classes                */
  int nx;      /* number of excluded frequency classes      */
  int *ax;     /* array of excluded frequency classes       */
  int E;       /* levels of exhaustive search               */
  int m;       /* maximum level searched exhaustively       */
  int x;       /* number of categories for cross-validation */
  int s;       /* seed for random number generator          */
  char **inputFiles;
  double u;    /* mutation rate */
  double c;    /* minimum change */
  int numInputFiles;
} Args;

typedef struct ints {
  int *a; /* array of integers */
  int  n; /* length of array   */
} Ints;

Args *getArgs(int argc, char *argv[]);
Args *newArgs();
void freeArgs(Args *args);
void printUsage();
void printSplash(char *version, char *date);

#endif
