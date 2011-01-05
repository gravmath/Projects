/*
  File: vector.cpp
  Author: David Fiske
  Purpose: General purpose vector class.

  $Header: /group/grt/project/cvsroot/newtonian_sph/vector.cpp,v 1.2 2001/12/05 23:50:11 drfiske Exp $
*/

#include "vector.h"

double length(const vector &v) {
  return sqrt(v.data[0]*v.data[0] + v.data[1]*v.data[1] + v.data[2]*v.data[2]);
}

void set(vector &v, int i, double x) {
  v.data[i-1] = x;
}
