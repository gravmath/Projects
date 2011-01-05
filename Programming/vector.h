/*
  File: vector.h
  Author: David Fiske
  Purpose: General purpose vector class.

  $Header: /group/grt/project/cvsroot/newtonian_sph/vector.h,v 1.2 2001/12/05 23:50:12 drfiske Exp $
*/

#ifndef __VECTOR__
#define __VECTOR__

#include <math.h>

class vector {
private:
  double data[3];
public:
  /**** Constructors ****/
  vector() {
    data[0] = 0.0;
    data[1] = 0.0;
    data[2] = 0.0;
  }
  vector(const vector &v) {
    data[0] = v.data[0];
    data[1] = v.data[1];
    data[2] = v.data[2];
  }
  vector(double x, double y, double z) {
    data[0] = x;
    data[1] = y;
    data[2] = z;
  }
  ~vector() {;}
  /**** Arithmatic ****/
  vector operator+(const vector &v) const {
    vector sum(data[0]+v.data[0], data[1]+v.data[1], data[2]+v.data[2]);
    return sum;
  }
  vector operator-(const vector &v) const {
    vector diff(data[0]-v.data[0], data[1]-v.data[1], data[2]-v.data[2]);
    return diff;
  }
  vector operator*(double alpha) const {
    vector prod(alpha*data[0],alpha*data[1],alpha*data[2]);
    return prod;
  }
  double operator*(const vector &v) const {
    return (data[0]*v.data[0] + data[1]*v.data[1] + data[2]*v.data[2]);
  }
  /**** Norms ****/
  friend double length(const vector&);
  /**** Assignment ****/
  vector& operator=(const vector &v) {
    data[0] = v.data[0];
    data[1] = v.data[1];
    data[2] = v.data[2];
    return *this;
  }
  vector& operator=(const double v[]) {
    data[0] = v[0];
    data[1] = v[1];
    data[2] = v[2];
    return *this;
  }
  vector& operator+=(const vector &v) {
    data[0] += v.data[0];
    data[1] += v.data[1];
    data[2] += v.data[2];
    return *this;
  }
  vector& operator-=(const vector &v) {
    data[0] -= v.data[0];
    data[1] -= v.data[1];
    data[2] -= v.data[2];
    return *this;
  }
  vector& operator*=(double alpha) {
    data[0] *= alpha;
    data[1] *= alpha;
    data[2] *= alpha;
    return *this;
  }
  vector& operator/=(double alpha) {
    data[0] /= alpha;
    data[1] /= alpha;
    data[2] /= alpha;
  }
  /**** Access ****/
  double operator()(int i) const {return data[i-1];}
  void set(int i, double x) {data[i-1] = x;}
  friend void set(vector&,int,double);
  void zero(void) {
    data[0] = 0.0;
    data[1] = 0.0;
    data[2] = 0.0;
  }
};

#endif

