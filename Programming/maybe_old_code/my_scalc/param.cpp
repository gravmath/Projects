// Param.cpp
// Implementation of simulation parameters class
// Keith Watt 1999

#include "param.hpp"    //<================= Changed to lower case 

// (KW) I've added a flag (dynCase) that determines whether the static or dynamic case should be run. 
//		if dynCase=0, it will run the static case, any other value run the dynamic case.
 
Param::Param()
{
	max.Resize(1,4);
	delta.Resize(1,4);
   smoothLength0.Resize(1,4);
	particleCenter0.Resize(1,4);
	particleVelocity0.Resize(1,4);
   smoothLength1.Resize(1,4);
	particleCenter1.Resize(1,4);
	particleVelocity1.Resize(1,4);
	max.Set(0.,0);
	max.Set(0.,1);
	max.Set(0.,2);
	max.Set(0.,3);
	delta.Set(0.,0);
	delta.Set(0.,1);
	delta.Set(0.,2);
	delta.Set(0.,3);
	diffCoeff = 0.;
	convGoal = 0.;
	omega   = 0.;
	boundCond = 0;
	printStep = 1;
	maxL = 0;
	maxN = 0;
	re   = 0.;
	dynCase = 0;  // <=================== Added flag to determine static or dynamic case; 0=static, non-zero=dynamic

	mass0 = 1.0;
	smoothLength0.Set(0.,0);
	smoothLength0.Set(0.25,1);
	smoothLength0.Set(0.25,2);
	smoothLength0.Set(0.25,3);
	particleCenter0.Set(0.,0);
	particleCenter0.Set(-0.5,1);
	particleCenter0.Set(0.,2);
	particleCenter0.Set(0.,3);
	particleVelocity0.Set(1.0,0);
	particleVelocity0.Set(0.,1);
	particleVelocity0.Set(0.,2);
	particleVelocity0.Set(0.,3);

	mass1 = 1.0;
	smoothLength1.Set(0.,0);
	smoothLength1.Set(0.25,1);
	smoothLength1.Set(0.25,2);
	smoothLength1.Set(0.25,3);
	particleCenter1.Set(0.,0);
	particleCenter1.Set(0.5,1);
	particleCenter1.Set(0.,2);
	particleCenter1.Set(0.,3);
	particleVelocity1.Set(1.0,0);
	particleVelocity1.Set(0.,1);
	particleVelocity1.Set(0.,2);
	particleVelocity1.Set(0.,3);

	exponent = 3.0;

	cout << "Default Param object created.\n";

}

Param::Param(double maxt, int maxx, int maxy, int maxz, double deltat,
		double deltax, double deltay, double deltaz, double diffC, double convG, double omegaIn,
		int bound, int printS, int mL, int mN, double reIn,
		double mass0In, double mass1In, double smoothLengthX0, double smoothLengthY0, double smoothLengthZ0,
		double smoothLengthX1, double smoothLengthY1, double smoothLengthZ1,
		double x0, double y0, double z0, double x1, double y1, double z1,
		double vx0, double vy0, double vz0, double vx1, double vy1, double vz1, double expIn,
      double x_min, double x_max, double y_min, double y_max, double z_min, double z_max, int dynCaseIn)  // <== Added flag to determine static or dynamic case; 0=static, non-zero=dynamic

{
	max.Resize(1,4);
	delta.Resize(1,4);
   smoothLength0.Resize(1,4);
	particleCenter0.Resize(1,4);
	particleVelocity0.Resize(1,4);
   smoothLength1.Resize(1,4);
	particleCenter1.Resize(1,4);
	particleVelocity1.Resize(1,4);
	max.Set(maxt,0);
	max.Set((double)maxx,1);
	max.Set((double)maxy,2);
	max.Set((double)maxz,3);
	delta.Set(deltat,0);
	delta.Set(deltax,1);
	delta.Set(deltay,2);
	delta.Set(deltaz,3);
	diffCoeff = diffC;
	convGoal = convG;
	omega = omegaIn;
	boundCond = bound;
	printStep = printS;
	maxL = mL;
	maxN = mN;
	re = reIn;
	dynCase = dynCaseIn; // <=================== Added flag to determine static or dynamic case; 0=static, non-zero=dynamic

	mass0 = mass0In;
	smoothLength0.Set(0.,0);
	smoothLength0.Set(smoothLengthX0,1);
	smoothLength0.Set(smoothLengthY0,2);
	smoothLength0.Set(smoothLengthZ0,3);
	particleCenter0.Set(0.,0);
	particleCenter0.Set(x0,1);
	particleCenter0.Set(y0,2);
	particleCenter0.Set(z0,3);
	particleVelocity0.Set(1.0,0);
	particleVelocity0.Set(vx0,1);
	particleVelocity0.Set(vy0,2);
	particleVelocity0.Set(vz0,3);

	mass1 = mass1In;
	smoothLength1.Set(0.,0);
	smoothLength1.Set(smoothLengthX1,1);
	smoothLength1.Set(smoothLengthY1,2);
	smoothLength1.Set(smoothLengthZ1,3);
	particleCenter1.Set(0.,0);
	particleCenter1.Set(x1,1);
	particleCenter1.Set(y1,2);
	particleCenter1.Set(z1,3);
	particleVelocity1.Set(1.0,0);
	particleVelocity1.Set(vx1,1);
	particleVelocity1.Set(vy1,2);
	particleVelocity1.Set(vz1,3);

   minimums.Set(x_min,0);
   minimums.Set(y_min,1);
   minimums.Set(z_min,2);
   maximums.Set(x_max,0);
   maximums.Set(y_max,1);
   maximums.Set(z_max,2);

	exponent = expIn;

	cout << "User-specified Param object created. \n";

}
