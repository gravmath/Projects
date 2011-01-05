/*
  File: math2.h
  Author: David Fiske
  Purpose: More math functions.

  $Header: /group/grt/project/cvsroot/newtonian_sph/math2.h,v 1.2 2001/12/05 22:25:55 drfiske Exp $
*/

#ifndef __MATH2__
#define __MATH2__

#include <math.h>

inline double pi(void) {
  return 3.1415926535897932384626433832795;
}

inline int sgn(int i) {
  return (i >= 0) ? 1 : -1;
}

inline int sgn(double x) {
  return (x >= 0.0) ? 1 : -1;
}

int closest_int(double);

#endif
