/***** search.c ***********************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Fri Nov 23 15:00:29 2018
 **************************************************/
#include "sfs.h"
#include "popSizes.h"
#include "interface.h"
#include "config.h"

void cpK(int *k1, int *k2, int m) {
  for(int i = 1; i <= m; i++)
    k2[i] = k1[i];
}

PopSizes *searchLevels(Sfs *sfs, Args *args) {
  int *kd;    /* array of levels */
  short ex;   /* exhaustive search? */

  ex = args->E;
  iniConfig(ex);
  PopSizes *ps = newPopSizes(sfs);
  newton(sfs, ps, args);
  l = ps->l;
  la = l;
  cpK(ps->k, ka, ps->m);
  /* iterate over the possible number of levels, n */
  for(int m = 2; m <= sfs->n; m++) {
    while((kd = nextConfig(ex)) != NULL) {
      ld = compPopSizes(kd, m, sfs, ps);
      if(ld > la) {
	la = ld;
	cpK(kd, ka, m);
      }
    }
    if(la <= l + 2) {
      compPopSizes(k, m - 1, sfs, ps);
      freeConfig(ex);
      return ps
    }
    cpK(ka, k, m);
    l = la;
  }
  compPopSizes(k, m - 1, sfs, ps);
  freeConfig(ex);
  return ps;
}
