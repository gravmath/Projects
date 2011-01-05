// Param: Class to hold all simulation parameters

#ifndef PARAM_HPP
#define PARAM_HPP

#include "tensor.hpp"

#define ZERO_FIXED 0
#define HARM_DECOMP 1

class Param {

// (KW) I've added a flag (dynCase) that determines whether the static or dynamic case should be run. 
//		if dynCase=0, it will run the static case, any other value run the dynamic case.
//		WARNING: THIS IS AN UNLOCKED CHANGE!

public:
	tensor	max; 			    	// Number of grid points (t,x,y,z)
	tensor	delta;				// Grid spacing (t,x,y,z)
	double	diffCoeff;			// Diffusion Coefficient for relaxation routine
	double	convGoal;			// Convergence goal for relaxation routine
	double	omega;				// Rotational velocity of frame
	int		boundCond;			// Boundary condition
	int 	   printStep;			// Print interval
	int		maxL;				   // Maximum value of l (spherical harmonics)
	int		maxN;				   // Maximum value of n (Misner harmonics)
	double   re;					// Radius of extraction

	int dynCase;					// 0=static case, non-zero=dynamic case  // <=== Added flag to determine static or dynamic case; 0=static, non-zero=dynamic

	double exponent;			   // Exponent for smoothing kernel

	double mass0;			      // Mass of particle 0
	tensor smoothLength0;		// Smoothing length (x,y,z)
	tensor particleCenter0;		// Center of particle (x,y,z)
	tensor particleVelocity0;	// Coordinate velocity of particle (vx, vy, vz)

	double mass1;			      // Mass of particle 0
	tensor smoothLength1;		// Smoothing length (x,y,z)
	tensor particleCenter1;		// Center of particle (x,y,z)
	tensor particleVelocity1;	// Coordinate velocity of particle (vx, vy, vz)

   array minimums;
   array maximums;

	Param();
	Param(double maxt, int maxx, int maxy, int maxz, double deltat,
		double deltax, double deltay, double deltaz, double diffC, double convG, double omegaIn,
		int bound, int printS, int mL, int mN, double reIn,
		double mass0In, double mass1In, double smoothLengthX0, double smoothLengthY0, double smoothLengthZ0,
		double smoothLengthX1, double smoothLengthY1, double smoothLengthZ1,
		double x0, double y0, double z0, double x1, double y1, double z1,
		double vx0, double vy0, double vz0, double vx1, double vy1, double vz1, double expIn,
      double x_min, double x_max, double y_min, double y_max, double z_min, double z_max, int dynCaseIn);  // <== Added flag to determine static or dynamic case; 0=static, non-zero=dynamic


	Param& operator<=(const Param &rval);

private:

};

inline Param &Param::operator<=(const Param &rval)
{
	max <= rval.max;
	delta <= rval.delta;
	diffCoeff = rval.diffCoeff;
	convGoal = rval.convGoal;
	omega = rval.omega;
	boundCond = rval.boundCond;
	printStep = rval.printStep;
	maxL = rval.maxL;
	maxN = rval.maxN;
	re = rval.re;

	dynCase = rval.dynCase;  // <=================== Added flag to determine static or dynamic case; 0=static, non-zero=dynamic

	mass0 = rval.mass0;
	smoothLength0 <= rval.smoothLength0;
	particleCenter0 <= rval.particleCenter0;
	particleVelocity0 <= rval.particleVelocity0;

	mass1 = rval.mass1;
	smoothLength1 <= rval.smoothLength1;
	particleCenter1 <= rval.particleCenter1;
	particleVelocity1 <= rval.particleVelocity1;

	exponent = rval.exponent;

   minimums <= rval.minimums;
   maximums <= rval.maximums;

	return *this;
}

#endif	// PARAM_HPP
