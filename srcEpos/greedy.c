/***** greedy.c ***********************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Fri Nov 23 16:18:38 2018
 **************************************************/
#include <stdlib.h>
#include "eprintf.h"

/* compInt compares two integers and returns their difference.
 * Called by qsort.
 */
static int compInt (const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
}

/* setupGreedy provides storage and preparation for the greedy search */
int *setupGreedy(int m, int n, int *start) {
  int *available = (int *)emalloc((n + 1) * sizeof(int));
  for(int i = 1; i <= n; i++)
    available[i] = 1;
  for(int i = 1; i < m; i++)
    available[start[i]] = 0;
  return available;
}

/* nextGreedy computes the next configuration using the greedy approach.
 * m: number of levels
 * n: sample size
 * start: start configuration of m-1 levels
 * setup: setup computation?
 */
int *nextGreedy(int m, int n, int *start, short setup) {
  static int *available, *k;
  
  if(setup) {
    available = setupGreedy(m, n, start);
    k = (int *)emalloc((m + 2) * sizeof(int));
    return NULL;
  }
  for(int i = 1; i < m; i++)
    k[i] = start[i];
  for(int i = 3; i <= n; i++) {
    if(available[i]) {
      k[m] = i;
      available[i] = 0;
      k[m + 1] = n + 1;
      qsort(k+1, m, sizeof(int), compInt);
      return k;
    }
  }
  free(available);
  free(k);

  return NULL;
}

void testGreedy() {
  int *k, *kd;
  int m = 2, n = 10;
  
  k = (int *)emalloc((n + 1) * sizeof(int));
  k[1] = 2;
  k[2] = 5;
  k[3] = n + 1;
  nextGreedy(m, n, k, 1);
  while((kd = nextGreedy(m, n, k, 0)) != NULL) {
    printf("%d", kd[1]);
    for(int i = 2; i <= m + 1; i++)
      printf(" %d", kd[i]);
    printf("\n");
  }
  free(k);
}
