/***** test.c *************************************
 * Description: Collect testing code.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Sun Nov 25 11:31:05 2018
 **************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "interface.h"

void testGreedy();
void testExhaustive();
void testNewton();

void test(Args *args) {
  printf("*** Test Greedy Search ***\n");
  testGreedy();
  printf("*** Test Exhaustive Search ***\n");
  testExhaustive();
  printf("*** Test Newton Computation ***\n");
  testNewton();
  free(args);
}
