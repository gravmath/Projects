/*
  File: grid.h
  Author: David Fiske
  Purpose: To hold a grid structure for the SPH problem.

  $Header: /group/grt/project/cvsroot/newtonian_sph/grid.h,v 1.1.1.1 2001/11/28 23:03:42 drfiske Exp $
*/

#ifndef __GRID__
#define __GRID__

#include "sph.h"
#include <math.h>
#include "math2.h"
#include "grid_point.h"
#include "cloud.h"
#include <stdio.h>

class grid_point;
class cloud;

class grid {
private:
  grid_point *point;
  int number;
  int total_number;
  double spacing;
  double smoothing;
  double lbnd[3];
  double ubnd[3];
  int valid_index(int i) const {
    return ( (i >= 0) && (i < total_number) );
  }
  int array2index(int i, int j, int k) const {
    return i + number*(j + number*k);
  }
public:
  /**** Constructors ****/
  grid();
  grid(const grid&);
  grid(char*);
  ~grid();
  /**** Physics ****/
  double rho(const vector&) const;
  double pressure(double) const;
  double pressure(const vector&) const;
  double W(double r) const {
    return W_spline(r,smoothing);
  }
  vector d1_W(const vector &r1, const vector &r2) const {
    return d1_W_spline(r1,r2,smoothing);
  }
  int search_range(void) const {
    return (int) ceil(2.0*smoothing/spacing);
  }
  /**** Operations ****/
  void reset(void);
  int closest_point(const vector&) const;
  int loop_subset(int,int) const;
  void contribute_mass2point(const vector &r, double m) {
    point[closest_point(r)].increase_mass(m);
  }
  void set_densities(void);
  void set_pressures(void);
  void set_densities_and_pressures(void);
  /**** Access ****/
  grid_point& operator()(int i) const {
    return point[i];
  }
};

#endif
