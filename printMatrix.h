/***** printMatrix.h ******************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Oct 30 10:27:55 2017
 **************************************************/
#ifndef PRINTMATRIX
#define PRINTMATRIX

#include <gsl/gsl_linalg.h>
#include "sfs.h"

void printMatrix(gsl_matrix *m, gsl_vector *b);

#endif
