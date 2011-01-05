/*
  File: grid.cpp
  Author: David Fiske
  Purpose: To hold a grid structure for the SPH problem.

  $Header: /group/grt/project/cvsroot/newtonian_sph/grid.cpp,v 1.2 2001/11/28 23:24:51 drfiske Exp $
*/

#include "grid.h"

grid::grid(const grid &g) {
  int i;
  number = g.number;
  total_number = g.total_number;
  point = new grid_point[total_number];
  for(i = 0; i < total_number; i++) {
    point[i] = g.point[i];
  }
}

/*
  The following constructor supports initializing the grid from a file.
  The format of the file is assumed to be as follows:

------Top of File------
<number of grid points in each direction (int)>
<spacing between points (double)>
<smoothing length  (double)>
<lowerbound x> <lowerbound y> <lowerbound z>
------End of File------
*/

grid::grid(char *name) {
  int i;
  int j;
  int k;
  double r[3];
  FILE *input;
  input = fopen(name, "r");
  fscanf(input, " %d\n", &number);
  fscanf(input, " %lf\n", &spacing);
  fscanf(input, " %lf\n", &smoothing);
  fscanf(input, " %lf %lf %lf\n", lbnd, lbnd+1, lbnd+2);
  fclose(input);
  total_number = number * number * number;
  point = new grid_point[total_number];
  ubnd[0] = lbnd[0] + spacing*number;
  ubnd[1] = lbnd[1] + spacing*number;
  ubnd[2] = lbnd[2] + spacing*number;
  r[0] = lbnd[0];
  r[1] = lbnd[1];
  r[2] = lbnd[2];
  for(k = 0; k < number; k++) {
    r[2] += k*spacing;
    for(j = 0; j < number; j++) {
      r[1] += j*spacing;
      for(i = 0; i < number; i++) {
	r[0] += i*spacing;
	point[array2index(i,j,k)].set_position(r);
      }
    }
  }
}

grid::~grid() {
  if (total_number > 0) {
    delete [] point;
  }
}

double
grid::rho(const vector &r) const {
  int i;
  double sum = 0.0;
  for(i = 0; i < total_number; i++) {
    sum += point[i].m() * W(length(r-point[i].r()),smoothing);
  }
}

double
grid::pressure(double rho) const {
  return rho;
}

double
grid::pressure(const vector &r) const {
  return rho(r);
}

void
grid::reset(void) {
  int i;
  for(i = 0; i < number; i++) {
    point[i].reset();
  }
}

int
grid::closest_point(const vector& r) const {
  int i;
  int j;
  int k;
  i = closest_int((r(1)-lbnd[0])/spacing);
  j = closest_int((r(2)-lbnd[1])/spacing);
  k = closest_int((r(3)-lbnd[2])/spacing);
  return array2index(i,j,k);
}

int
grid::loop_subset(int center, int N) const {
  static int result = -1;
  static int i;
  static int j;
  static int k;
  if (result == -1) {
    i = -N;
    j = -N;
    k = -N;
  }
  else {
    if (i < N) {
      i++;
      result = center + array2index(i,j,k);
      result = (valid_index(result)) ? result : loop_subset(center,N);
    }
    else {
      i = -N;
      if (j < N) {
	j++;
	result = center + array2index(i,j,k);
	result = (valid_index(result)) ? result : loop_subset(center,N);
      }
      else {
	j = -N;
	if (k < N) {
	  k++;
	  result = center + array2index(i,j,k);
	  result = (valid_index(result)) ? result : loop_subset(center,N);
	}
	else {
	  result = -1;
	}
      }
    }
  }
  return result;
}

void
grid::set_densities(void) {
  int i;
  for(i = 0; i < number; i++) {
    point[i].set_density(rho(point[i].r()));
  }
}

void
grid::set_pressures(void) {
  int i;
  for(i = 0; i < number; i++) {
    point[i].set_pressure(pressure(point[i].rho()));
  }
}

void
grid::set_densities_and_pressures(void) {
  int i;
  double rho_temp;
  double pressure_temp;
  for(i = 0; i < number; i++) {
    rho_temp = rho(point[i].r());
    pressure_temp = pressure(rho_temp);
    point[i].set_density(rho_temp);
    point[i].set_pressure(pressure_temp);
    point[i].set_ratio();
  }
}
