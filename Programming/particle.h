/*
  File: particle.h
  Author: David Fiske
  Purpose: To represent a particle for SPH.

  $Header: /group/grt/project/cvsroot/newtonian_sph/particle.h,v 1.1.1.1 2001/11/28 23:03:42 drfiske Exp $
*/

#ifndef __PARTICLE__
#define __PARTICLE__

#include "vector.h"
#include "grid.h"

class vector;
class grid;

class particle {
private:
  double mass;
  vector position;
  vector velocity;
  vector acceleration;
public:
  /**** Constructors ****/
  particle();
  particle(const particle &p) {
    mass = p.mass;
    position = p.position;
    velocity = p.velocity;
  }
  particle(double m, const vector &r, const vector &v) {
    mass = m;
    position = r;
    velocity = v;
  }
  ~particle() {;}
  /**** Operations ****/
  void move(const grid&);
  particle& operator=(const particle &p) {
    mass = p.mass;
    position = p.position;
    velocity = p.velocity;
    return *this;
  }
  /**** Data Access ****/
  void set_data(double m, double r[3], double v[3]) {
    mass = m;
    position = r;
    velocity = v;
  }
  double m(void) {
    return mass;
  }
  vector& r(void) {
    return position;
  }
  double r(int i) {
    return position(i);
  }
  vector& v(void) {
    return velocity;
  }
  double v(int i) {
    return velocity(i);
  }
};

#endif
