Param

1) Added line to properly resize the particle arrays

   smoothLength0.Resize(1,4);
	particleCenter0.Resize(1,4);
	particleVelocity0.Resize(1,4);
   smoothLength1.Resize(1,4);
	particleCenter1.Resize(1,4);
	particleVelocity1.Resize(1,4);

MSfield

1) Added private member data pointer to Particle - default to NULL
2)	Added public member function RegisterParticles
3)	Changed all calls to internal matterVar
4) Rewrote calcDgDphi so that it returns a tensor that can be contracted
   with the result of Particle::matterVar
5) added private tensor member data  matter_piece and field_piece
6) added the following lines to the MSField constructor
	   my_particle_partner = NULL;
	   matter_piece.Resize(2,4,4);
	   field_piece. Resize(2,4,4);
	   g.           Resize(2,4,4);
7) added dgdphi tensor to the private member data and the following to ctor
      dgdphi.Resize(2,4,4);
8) added member data omega to the private data to facilitate speed rather than accessing
   simParam.omega each time



Particle

1)	Added private member data pointer to MSField - default to NULL
2) Added public member function RegisterField

Param
1)	Added parameters to specify the minimum and maximum extent of the grid in
   each of the three directions


