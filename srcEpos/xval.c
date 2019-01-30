/***** xval.c *************************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Wed Dec 19 17:17:04 2018
 **************************************************/
#include <time.h>
#include "sfs.h"
#include "eprintf.h"
#include "util.h"

SfsSet *newSfsSet(Sfs *sfs, Args *args) {
  int x = args->x;
  int n = sfs->n;
  SfsSet *ss = (SfsSet *)emalloc(sizeof(SfsSet));
  ss->n = x;
  ss->train = (Sfs **)emalloc(x * sizeof(SfsSet *));
  ss->test  = (Sfs **)emalloc(x * sizeof(SfsSet *));
  for(int i = 0; i < x; i++) {
    Sfs *s = newSfs(n, args);
    s->a = sfs->a;
    s->l = sfs->l;
    ss->train[i] = s;
    s = newSfs(n, args);
    s->a = sfs->a;
    s->l = sfs->l;
    ss->test[i]  = s;
  }
  return ss;
}

void freeSfsSet(SfsSet *ss) {
  for(int i = 0; i < ss->n; i++) {
    freeSfs(ss->train[i]);
    freeSfs(ss->test[i]);
  }
  free(ss->train);
  free(ss->test);
  free(ss);
}

/* addSfs adds the SFS of b to that of a */
void addSfs(Sfs *a, Sfs *b) {
  for(int i = 0; i <= a->a; i++) {
    a->G[i] += b->G[i];
  }
  a->p += b->p;
  a->l += b->l;
}

/* splitSfs splits the site frequency spectrum sfs into x parts
 * ready for cross-validation.
 */
SfsSet *splitSfs(Sfs *sfs, Args *args, gsl_rng *r) {
  SfsSet *ss = newSfsSet(sfs, args);
  if(args->x == 1) { /* no cross-validation */
    freeSfs(ss->train[0]);
    freeSfs(ss->test[0]);
    ss->train[0] = copySfs(sfs);
    ss->test[0] = copySfs(sfs);
    return ss;
  }
  int *a = (int *)emalloc(sfs->p * sizeof(int));
  long n = 0;
  for(int i = 1; i <= sfs->a; i++) {
    for(int j = 0; j < sfs->G[i]; j++) {
      a[n++] = i;
    }
  }
  shuffle(a, n, r);
  n /= args->x;
  long x = 0;
  long l = sfs->l / args->x;
  for(int i = 0; i < args->x; i++) {
    Sfs *ns = ss->test[i];
    for(long j = 0; j < n; j++) {
      ns->G[a[x]]++;
      ns->p++;
      x++;
    }
    ns->G[0] = l - ns->p;
    ns->l = l;
  }
  for(int i = 0; i < args->x; i++) {
    freeSfs(ss->train[i]);
    ss->train[i] = copySfs(ss->test[i]);
    for(int j = 0; j < args->x; j++) {
      if(i != j)
	addSfs(ss->train[i], ss->test[j]);
    }
  }
  for(int i = 0; i < args->x; i++) {
    for(int j = 0; j < args->nx; j++) {
      ss->train[i]->G[args->ax[j]] = -1;
      ss->test[i]->G[args->ax[j]]  = -1;
    }
  }
  free(a);

  return ss;
}
