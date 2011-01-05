// MSField.h
// Class which encapsulates a gridded Misner scalar field

#ifndef MSFIELD_HPP
#define MSFIELD_HPP

//#define PART			<===================== Uncomment to use Particle code
#include <iostream.h>
#include "tensor.hpp"
#include "param.hpp"
#include "constant.hpp"
#include "util.hpp"

#ifdef PART
#include "particle.hpp"
class Particle;
#endif

class MSField {

public:
	// Create a flat-space field
	MSField();

	//Destructor
	~MSField();

// Retrieve the value of harmonic coefficient for the requested harmonic
	double HarmCoeff(int c);
	// Retrieve the entire harmonic coefficient array
	array HarmCoeff();

// Retrieve the value of Phi at time step q and grid position i,j,k
	double Phi(int q, int i, int j, int k);
	// Retrieve the value of Pi at time step q and grid position i,j,k
	double Pi(int q, int i, int j, int k);

	// Retrieve simulation parameters
    Param Params();

	// Set the value of Phi at time step q and grid position i,j,k
	void SetPhi(double value, int q, int i, int j, int k);
	// Set the value of Pi at time step q and grid position i,j,k
	void SetPi(double value, int q, int i, int j, int k);

	// Set the simulation parameters
	void SetParams(Param simParamsIn);

	// Calculate the static-source omega
	double calcOmega(double currOmega);
    // Calculate the force applied by the motor
	double calcMotorForce(int i);

	// Calculate the Phi time derivative
	double dPhidt(double stepFrac, int q, int i, int j, int k);
	// Calculate the Pi time derivative
	double dPidt(double stepFrac, double d, int q, int i, int j, int k);

	// (KW) Added method to calculate the static source term; removed other two methods (matterVar and calcDensity)
	double calcStaticSource(int a, int x, int y, int z);

	// Calculate the bare mass
	double calcBareMass();

	// Calculate the metric tensor
	tensor GetInvG(int q, int i, int j, int k);
	double metricG(int q, int i, int j, int k, int a, int b);

	// Calculate dg/dphi
	double calcDgDphi(int q, int i, int j, int k,  int a, int b);
    tensor calcDgDphi(int q, int i, int j, int k);   

	// Calculate the field energy
	// q is the time step for which the energy is to be computed
	double calcEnergy(int q);

	// Calculate the dipole portion of the Poynting flux
	double calcDipolePoyntingFlux();

	// Calculate the quadrupole portion of the Poynting flux
	double calcQuadPoyntingFlux();

	// Initialize the field with flat-space values
	void flatSpace(void);

	// Relax the grid to achieve an initial solution, returns the computed value of omega
	double relax(void);

	// Perform a 2nd-order Runge-Kutta step
	// stepFrac is the percentage of the system timestep, simParams.delta.Val(0)
	void rk2(double stepFrac, double d);

	// Calculate the boundary values based on the boundary conditions
	// q is the time step for which the boundary is to be computed
	void calcBoundary(int q);

	// Calculate the position on the grid based on grid indices
	tensor Pos(int i, int j, int k);
	double Pos(int comp, int i, int j, int k);
	double rDist(int i, int j, int k);

	// Calculate the nearest grid index based on a position
	int Index(int component, double x);


	// These functions are used for the harmonic decomposition boundary conditions
	
	// Calculate l, m, and n from a harmonic index
	tensor getK(int k);
	int getK(int comp, int k);
	// Calculate a harmonic index from l, m, and n
	int getK(int l, int m, int n);
	// Calculate the total number of harmonic indices
	int numHarm(void);
	// Calculate the normalization tensor
	tensor calcNormTensor(void);
	// Calculate the harmonic decomposition
	void calcMisnerHarmonics(void);
	//Calculate the field from a harmonic decomposition
	double restoreFromMisnerHarmonics(tensor pos);
	double restoreFromMisnerHarmonics(int i, int j, int k);
	// Initialize the harmonic decomposition
	void initMisnerHarmonics(void);

#ifdef PART
	void RegisterParticles(Particle *particle);
#endif

private:
	double ****phi;  //declare phi as a quad pointer to doubles
    double ****pi;  //declare pi as a quad pointer to doubles
	double *b;
	tensor g;
	tensor dgdphi;        //added by CS 2/5/00
	tensor g_inv;         //added by CS 2/9/00
	tensor four_pos;      //added by CS 2/9/00
	double omega;
	Param simParams;

#ifdef PART
	Particle *my_particle_partner;
	tensor matter_piece;  //added by CS 2/5/00
	tensor field_piece;   //added by CS 2/5/00
#endif
};

// The destructor needs to return the memory used for the
//	dynamically-allocated Particle array
inline MSField::~MSField()
{

}

// Return the value of the harmonic coefficientfor the specified harmonic
//	This routine is only needed for exterior classes to access the
//	private member data.
inline double MSField::HarmCoeff(int c)
{
	return b[c];
}

// Return the value of the field at the requested grid location
//	and for the requested timeslice (0="old", 1="current", 2="new")
//	This routine is only needed for exterior classes to access the
//	private member data.
inline double MSField::Phi(int t, int x, int y, int z)
{
	return phi[t][x][y][z];
}

// Return the value of the field time derivative at the requested grid 
//	location and for the requested timeslice (0="old", 1="current", 2="new")
//	This routine is only needed for exterior classes to access the
//	private member data.
inline double MSField::Pi(int t, int x, int y, int z)
{
	return pi[t][x][y][z];
}

// Set the value of the field at the requested grid location
//	and for the requested timeslice (0="old", 1="current", 2="new")
//	This routine is only needed for exterior classes to access the
//	private member data.
inline void MSField::SetPhi(double value, int t, int x, int y, int z)
{
	phi[t][x][y][z] = value;
	return;
}

// Set the value of the field time derivative at the requested grid 
//	location and for the requested timeslice (0="old", 1="current", 2="new")
//	This routine is only needed for exterior classes to access the
//	private member data.
inline void MSField::SetPi(double value, int t, int x, int y, int z)
{
	pi[t][x][y][z] = value;
	return;
}


// Return the simulation paramters.
//	This routine is only needed for exterior classes to access the
//	private member data.
inline Param MSField::Params(void)
{
	return simParams;
}

// Return a single component of the coordinate position of a grid-indexed location
// This was done to avoid having to create a tensor object.  In actual usage, you rarely 
// actually need the position as a tensor anyway - you almost always want the coordinates.
inline double MSField::Pos(int comp, int i, int j, int k)
{
	if (comp==1)
	{
		return ((double)(i)-0.5*(simParams.max.Val(1)-1.0))*simParams.delta.Val(1);
	}
	else if (comp==2)
	{
		return ((double)(j)-0.5*(simParams.max.Val(2)-1.0))*simParams.delta.Val(2);
	}
	else if (comp==3)
	{
		return ((double)(k)-0.5*(simParams.max.Val(3)-1.0))*simParams.delta.Val(3);
	}
	else
	{
		cout << "ERROR: Cannot determine position for component " << comp << "\n";
		return 0.;
	}
}

inline double MSField::rDist(int i, int j, int k)
{
	return sqrt(Pos(1,i,j,k)*Pos(1,i,j,k) + Pos(2,i,j,k)*Pos(2,i,j,k) + Pos(3,i,j,k)*Pos(3,i,j,k));
}


#endif	// MSFIELD_HPP



