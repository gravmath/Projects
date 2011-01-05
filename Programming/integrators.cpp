/*
  File: integrators.cpp
  Author: David Fiske
  Purpose: General purpose integration routines.

  $Header: /group/grt/project/cvsroot/newtonian_sph/integrators.cpp,v 1.1.1.1 2001/11/28 23:03:42 drfiske Exp $
*/

#include "integrators.h"

double euler_integrate(double y0, double dy_dt, double dt) {
  return y0 + dt*dy_dt;
}
