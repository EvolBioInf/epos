/***** foldedE.h **********************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Feb 26 16:30:12 2018
 **************************************************/
#ifndef FOLDED
#define FOLDED

#include "popSizes.h"
#include "sfs.h"
#include "interface.h"
#include <gsl/gsl_linalg.h>

int foldedEcompPopSizes(Sfs *sfs, PopSizes *ps, Args *args);
double foldedEpsi(PopSizes *ps, Sfs *sfs);
gsl_vector *foldedEgetResVect(Sfs *sfs, PopSizes *ps);
gsl_matrix *foldedEgetCoeffMat(Sfs *sfs, PopSizes *ps, Args *args);

#endif
