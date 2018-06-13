/***** printMatrix.c ******************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Oct 30 10:18:44 2017
 **************************************************/
#include <gsl/gsl_linalg.h>
#include <stdio.h>
#include "sfs.h"

void printMatrix(gsl_matrix *m, gsl_vector *b){
  int i, j;
  printf("%s\n\n", "Data in LAPAC-Example format");
  printf("%d %d\n\n", (int)m->size1, 1);
  for(i=0;i<m->size1;i++){
    printf("%.3e", gsl_matrix_get(m,i,0));
    for(j=1;j<m->size2;j++)
      printf(" %.3e", gsl_matrix_get(m,i,j));
    printf("\n");
  }
  printf("\n");
  if(!b)
    return;
  for(i=0;i<b->size;i++)
    printf("%.3e\n",gsl_vector_get(b,i));
}

