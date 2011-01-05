/*
  File: sph.cpp
  Author: David Fiske
  Purpose: To implement sph-specific mathematical functions.

  $Header: /group/grt/project/cvsroot/newtonian_sph/sph.cpp,v 1.1.1.1 2001/11/28 23:03:42 drfiske Exp $
*/

#include "sph.h"

double W_spline(double r, double h) {
  double v = r/h;
  double c = pi()*h*h*h;
  if (v <= 1.0) {
    return (1.0 + v*v * (-1.5 + 0.75*v))/c;
  }
  else if (v <= 2.0) {
    return 0.25*(2.0 - v)*(2.0 - v)*(2.0 - v)/c;
  }
  else {
    return 0.0;
  }
}

vector d1_W_spline(const vector &r1, const vector &r2, double h) {
  vector c = (r1-r2);
  double r = length(c);
  double v = r/h;
  c *= 2.0/(pi()*h*h*h*h*r);
  if (v <= 1.0) {
    c *= (3.0*v)*(1 + 0.75*v);
    return c;
  }
  else if (v <= 2.0) {
    c *= -0.75*(2.0 - v)*(2.0 - v);
    return c;
  }
  else {
    c *= 0.0;
    return c;
  }
}
