/***** exhaustive.c *******************************
 * Description: Knuth's Algorithm L for finding all
 * combinations; Vol 4a of TACP.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Fri Nov 23 16:18:33 2018
 **************************************************/
#include <stdlib.h>
#include "eprintf.h"

/* se sets up the computation of all n \choose m combinations */
static int *se(int m, int n) {
  int *c;

  c = (int *)emalloc((m + 3) * sizeof(int));
  for(int i = 1; i <= m; i++)
    c[i] = i - 1;
  c[m + 1] = n;
  c[m + 2] = 0;
  
  return c;
}

/* ne returns the next combination in a sequence of of n \choose m combinations */
static int *ne(int m, int n, int *k, short setup) {
  static short finished = 0;
  static int *c;

  if(setup) {
    finished = 0;
    c = se(m, n);
    return NULL;
  }
  if(finished) {
    free(c);
    return NULL;
  }
  k[1] = 2;
  for(int i = 1; i <= m; i++)
    k[i + 1] = c[i] + 3;
  k[m + 2] = n + 3;
  int j = 1;
  while(c[j] + 1 == c[j + 1]) {
    c[j] = j - 1;
    j++;
  }
  if(j > m) {
    finished = 1;
  }
  if(!finished)
    c[j]++;

  return k;
}

int *nextExhaustive(int m, int n, int *k, short setup) {
  return ne(m-1, n - 2, k, setup);
}

void testExhaustive() {
  int n = 5, m = 3;
  int *k = (int *)emalloc((m + 2) * sizeof(int));

  nextExhaustive(m, n, k, 1);
  while((k = nextExhaustive(m, n, k, 0)) != NULL) {
    printf("%d", k[1]);
    for(int i = 2; i <= m + 1; i++)
      printf(" %d", k[i]);
    printf("\n");
  }
}


