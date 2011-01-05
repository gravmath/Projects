/*
  File: sph.h
  Author: David Fiske
  Purpose: To implement sph-specific mathematical functions.

  $Header: /group/grt/project/cvsroot/newtonian_sph/sph.h,v 1.1.1.1 2001/11/28 23:03:42 drfiske Exp $
*/

#ifndef __SPH__
#define __SPH__

#include "math2.h"
#include "vector.h"

/**** The smoothing kernal ****/
double W_spline(double,double);
vector d1_W_spline(const vector&,const vector&,double);

#endif
