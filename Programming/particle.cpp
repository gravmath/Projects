/*
  File: particle.cpp
  Author: David Fiske
  Purpose: To represent a particle for SPH.

  $Header: /group/grt/project/cvsroot/newtonian_sph/particle.cpp,v 1.2 2001/11/30 17:53:24 drfiske Exp $
*/

#include "particle.h"

#define TIME_SPACING 0.1

void
particle::move(const grid &g) {
  int i;
  int center = g.closest_point(position);
  int search = g.search_range();
  double my_rho = g.rho(position);
  double my_pressure = g.pressure(my_rho);
  double ratio = my_pressure/(my_rho*my_rho);
  /*
    dv/dt = - Sum_j m_j (P_i/(rho_i^2) + P_j/(rho_j^2)) grad W
            + m_i (P_i/(rho_i^2) + P_closest/(rho_closest^2)) grad W
  */
  acceleration = g.d1_W(position,g(center).r()) 
    * mass * (ratio + g(center).ratio());
  for(i = g.loop_subset(center,search); i != -1; 
      i = g.loop_subset(center,search)) {
    acceleration -= g.d1_W(position,g(i).r()) 
      * g(i).m() * (ratio + g(i).ratio());
  }
  for(i = 1; i <= 3; i++) {
    position.set(i,euler_integrate(position(i),velocity(i),TIME_SPACING));
  }
  for(i = 1; i <= 3; i++) {
    velocity.set(i,euler_integrate(velocity(i),acceleration(i),
				      TIME_SPACING));
  }
}
