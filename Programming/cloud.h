/*
  File: cloud.h
  Author: David Fiske
  Purpose: To implement matter as an array of particles with an SPH
  approximation.

  $Header: /group/grt/project/cvsroot/newtonian_sph/cloud.h,v 1.1.1.1 2001/11/28 23:03:42 drfiske Exp $
*/

#ifndef __CLOUD__
#define __CLOUD__

#include "vector.h"
#include "particle.h"
#include "grid.h"
#include <stdio.h>

class vector;
class particle;
class grid;

class cloud {
private:
  particle *matter;
  int number;
public:
  /**** Constructors ****/
  cloud() {;}
  cloud(const cloud&);
  cloud(int);
  cloud(char*);
  ~cloud();
  /**** Operations ****/
  void transfer_mass2grid(grid&);
  /**** Physics ****/
  void evolve(const grid&);
};

#endif
