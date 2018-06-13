/***** xval.h *************************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Mar 12 16:08:30 2018
 **************************************************/
#ifndef XVAL_H
#define XVAL_H

#include "sfs.h"
#include "interface.h"
#include "gsl_rng.h"

void xvalM(Sfs *sfs, Args *args, gsl_rng *r);

#endif
