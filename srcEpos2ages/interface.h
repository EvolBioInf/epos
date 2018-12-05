#ifndef INTERFACE
#define INTERFACE

/* define argument container */
typedef struct args{
  char h;   /* help message?    */
  char v;   /* version message? */
  char e;   /* error message?   */
  int  n;   /* sample size      */
  char **inputFiles;
  int numInputFiles;
} Args;

Args *getArgs(int argc, char *argv[]);
void printUsage();
void printSplash(char *version, char *date);

#endif
