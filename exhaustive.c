/***** exhaustive.c *******************************
 * Description: Generate all combinations of 
 * levels using Knuth's Algorithm T.
 * Reference: Donald E. Knuth (2011). The Art of
 *   Computer Programming. Volume 4A: Combinatorial
 *   Algorithms, Part 1. Addison Wesley, p. 359.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Fri Nov 23 16:18:33 2018
 **************************************************/
#include <stdlib.h>
#include "eprintf.h"

/* se sets up the computation of all n \choose m combinations */
static int *se(int m, int n) {
  int *c;
  /* T1. [Initialize.] */
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
  static int j;
  int x;

  if(setup) {
    finished = 0;
    c = se(m, n);
    j = m;
    return NULL;
  }
  if(finished) {
    free(c);
    return NULL;
  }
  /* T2. [Visit.] */
  k[1] = 2;
  for(int i = 1; i <= m; i++)
    k[i + 1] = c[i] + 3;
  k[m + 2] = n + 3;
  if(j > 0) {
    x = j;
  } else {
    /* T3. [Easy case?] */
    if(c[1] + 1 < c[2]) {
      c[1]++;
      return k;
    } else {
      j = 2;
    }
    /* T4. [Find j.] */
    c[j - 1] = j - 2;
    x = c[j] + 1;
    while(x == c[j + 1]) {
      j++;
      c[j - 1] = j - 2;
      x = c[j] + 1;
    }
    /* T5. [Done?] */
    if(j > m) {
      finished = 1;
    }
  }
  /* T6. [Increase c_j.] */
  c[j] = x;
  j--;

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


