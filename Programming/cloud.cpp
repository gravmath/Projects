/*
  File: cloud.cpp
  Author: David Fiske
  Purpose: To implement matter as an array of particles with an SPH
  approximation.

  $Header: /group/grt/project/cvsroot/newtonian_sph/cloud.cpp,v 1.2 2001/11/28 23:24:50 drfiske Exp $
*/


#include "cloud.h"

cloud::cloud(const cloud &c) {
  int i;
  number = c.number;
  matter = new particle[number];
  for(i = 0; i < number; i++) {
    matter[i] = c.matter[i];
  }
}

cloud::cloud(int n) {
  number = n;
  matter = new particle[number];
}

/*
  The following constructor supporst initializing the matter cloud from a file.
  The format of the file is assumed to be as follows.

-------Top of File------
<number of particles (int)>
<mass of particle 1 (double)> (<x> <y> <z>) (<vx> <vy> <vz>)
<mass of particle 2 (double)> (<x> <y> <z>) (<vx> <vy> <vz>)
...
<mass of particle N (double)> (<x> <y> <z>) (<vx> <vy> <vz>)
-------End of File------
*/

cloud::cloud(char *name) {
  FILE *input;
  int i;
  double m;
  double r[3];
  double v[3];
  input = fopen(name, "r");
  fscanf(input, " %d\n", &number);
  matter = new particle[number];
  for(i = 0; i < number; i++) {
    fscanf(input, " %lf ( %lf %lf %lf ) ( %lf %lf %lf )\n", &m,
	   r, r+1, r+2, v, v+1, v+2);
    matter[i].set_data(m,r,v);
  }
  fclose(input);
}

cloud::~cloud() {
  if (number != 0) {
    delete [] matter;
  }
}

void
cloud::transfer_mass2grid(grid &g) {
  int i;
  for(i = 0; i < number; i++) {
    g.contribute_mass2point(matter[i].r(), matter[i].m());
  }
}

void
cloud::evolve(const grid &g) {
  int i;
  for(i = 0; i < number; i++){
    matter[i].move(g);
  }
}
