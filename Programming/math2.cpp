/*
  File: math2.cpp
  Author: David Fiske
  Purpose: More math functions.

  $Header: /group/grt/project/cvsroot/newtonian_sph/math2.cpp,v 1.2 2001/12/05 22:25:55 drfiske Exp $
*/

#include "math2.h"

int closest_int(double x) {
  double integer_part;
  double fractional_part;
  int i;
  fractional_part = modf(x, &integer_part);
  i = (int) integer_part;
  if (fabs(fractional_part) >= 0.5) {
    if (i > 0) {
      i++;
    }
    else if (i < 0) {
      i--;
    }
    else {
      i = sgn(fractional_part);
    }
  }
  return i;
}

