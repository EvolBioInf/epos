/***** unfolded.h *********************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Feb 26 16:33:19 2018
 **************************************************/
#ifndef UNFOLDEDH
#define UNFOLDEDH

#include "sfs.h"
#include "popSizes.h"
#include "interface.h"
#include <gsl/gsl_linalg.h>

void unfoldedCompPopSizes(Sfs *sfs, PopSizes *ps, Args *args);
double unfoldedPsi(PopSizes *ps, Sfs *sfs);
gsl_vector *unfoldedGetResVect(Sfs *sfs, PopSizes *ps);
gsl_matrix *unfoldedGetCoeffMat(Sfs *sfs, PopSizes *ps, Args *args);

#endif
