/*
  File: grid_point.h
  Author: David Fiske
  Purpose: To define a class containing the information required by a 
  grid point.

  $Header: /group/grt/project/cvsroot/newtonian_sph/grid_point.h,v 1.1.1.1 2001/11/28 23:03:42 drfiske Exp $
*/

#ifndef __GRID_POINT__
#define __GRID_POINT__

#include "vector.h"
#include "particle.h"
#include "cloud.h"

class vector;
class particle;
class cloud;

class grid_point {
private:
  vector position;
  double density;
  double mass;
  double pressure;
  double p_rho_ratio;
public:
  /**** Constructors ****/
  grid_point() {;}
  grid_point(const grid_point &p) {
    position = p.position;
    density = p.density;
    mass = p.mass;
    pressure = p.pressure;
  }
  grid_point(const vector &v) {
    position = v;
  }
  ~grid_point() {;}
  /**** Operations ****/
  void reset(void) {
    mass = 0.0;
  }
  void increase_mass(double m) {
    mass += m;
  }
  void set_density(double rho) {
    density = rho;
  }
  void set_pressure(double p) {
    pressure = p;
  }
  void set_ratio(void) {
    p_rho_ratio = pressure/(density*density);
  }
  /**** Data Access ****/
  void set_position(double r[3]) {
    position = r;
  }
  vector& r(void) {
    return position;
  }
  double m(void) {
    return mass;
  }
  double rho(void) {
    return density;
  }
  double p(void) {
    return pressure;
  }
  double ratio(void) {
    return p_rho_ratio;
  }
};

#endif
