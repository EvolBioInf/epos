/***** config.c ***********************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Fri Nov 23 16:03:40 2018
 **************************************************/
#include "exhaustive.h"
#include "greedy.h"


void iniConfig(short ex, int n, int k) {
  if(ex)
    iniExhaustive(n, k);
  else
    iniGreedy(n, k);
}

int *nextConfig(short ex) {
  if(ex)
    return nextExhaustive();
  else
    return nextGreedy();
}
