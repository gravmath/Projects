/*
  File: init.hpp
  Author: Conrad Schiff
  Purpose: Routine to initialize N particles for Newtonian SPH.
*/

#ifndef __INIT__
#define __INIT__

#include <stdlib.h>


class init {
private:
  indouble mass;
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