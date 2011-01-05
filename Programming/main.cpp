/*
  File: main.cpp
  Author: David Fiske
  Purpose: A sample main for Newtonian SPH.

  $Header: /group/grt/project/cvsroot/newtonian_sph/main.cpp,v 1.1.1.1 2001/11/28 23:03:43 drfiske Exp $
*/

#include "grid.h"
#include "cloud.h"

int main(void) {
  grid g("grid1.sph");
  cloud star("star1.sph");
  star.transfer_mass2grid(g);
  g.set_densities_and_pressures();
  star.evolve(g);
  return 1;
}
