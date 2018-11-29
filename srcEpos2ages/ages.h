/***** ages.h *************************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Thu Nov 29 12:31:26 2018
 **************************************************/
#ifndef AGES
#define AGES

#include "popSizes.h"

typedef struct ages {
  int n;     /* sample size */
  double *N; /* population sizes N[2], N[3],..., N[n] */
  double *a; /* ages        */
  double *s; /* sizes       */
} Ages;

Ages *compAges(PopSizes *ps);
void freeAges(Ages *a);
void printAges(Ages *a);

#endif
