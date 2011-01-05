/* Main program driver */
// This is the changed versin

#include "msfield.hpp"   //<====================== Changed to lower case

int main(int argc, char *argv[])
{
	int i,j,k;

	// Load the simulation parameters from a file (or hard code here)

	// (KW) Added flag to Param object which determines if the static or dynamic case is being evolved.

	// Set the grid spacing.  Timestep is 1/4 grid spacing. Units are r/M.
	double gridspacing = 0.5;

	// Set the starting orbital seperation.  NS are on x-axis. Units are r/M.
	double R = 7.5;

	// Calculate the number of gridpoints between the center of the particle and the center of the grid
	double nR = R/gridspacing;

	// Set the radius of the NS. The constant is the number of gridpoints in each direction.  Units are r/M.
	double radius = 0.3;

	// Set the central masses of one NS.
	double mass0 = 0.1;

	// Set the mass ratio between the NS
	double mu = 1.;

	// Set the exponent of the smoothing kernel
	double exponent = 3.0;

	// Set the number of gridpoints
	double gridpoints = 39.; // Increase 2*(R+radius)/gridspacing + 6 to the next odd integer

	// Set the radius of extraction
	double re = 9.; // Increase R+radius+(1.5*gridspacing) to next integer gridpoint;

	// Calculate the number of gridpoints between the radius of extraction and the center of the grid
	double nRe = re/gridspacing;

	// Set the ending orbital separation.  By default, this is when the NS just touch. Units are r/M.
	double Rend = radius;

	// Set the number of snapshots to take of the evolution
	int snapshots = 10;

	// Calculate the ending grid spacing.  Units are r/M.
	double endgridspacing = Rend/nR;

	// Calculate the delta gridspacing
	double deltaspacing = (gridspacing + endgridspacing)/((double)snapshots);

	// Set initial omega
	double startOmega = 0.1;

	// Set convergence tolerance
	double convCriterion = 0.001;

	// Set diffusion coefficient
	double diffCoeff = 2.42;

	// Calculate the starting size of the grid.  Radius of extraction is R+radius+3*gridspacing3
	// Create a simulation parameters object
	Param simParams(0., (int)gridpoints, (int)gridpoints, (int)gridpoints, gridspacing/4., gridspacing, gridspacing, gridspacing,
			diffCoeff, convCriterion, startOmega, HARM_DECOMP, 500, 2, 2, re,
			mass0, mu*mass0, radius,radius,radius, radius,radius,radius,
			-R,0.0,0.0, R,0.0,0.0, 0.,0.,0., 0.,0.,0.,exponent,
         -1.3,-1.3,-1.3,1.3,1.3,1.3,0);   // <======================= Added 0 as last parameter, set to non-zero for dynamic case
										  // <======================= According to my reading of the simParams() routine, this should
	                                      //                            be -1.3,1.3,-1.3,1.3,-1.3,1.3,0);

	// Create a flat scalar field and send it a copy of the simulation parameters
	MSField field;
    field.SetParams(simParams);
 	cout << "Field object created.\n";

#ifdef PART
	Particle swarm;
	swarm.SetParams(simParams);  <======================== Seems to crash when attempting to set inverse deltas

	cout << "Particle params initialized.\n";

   //Register the field and particles with each other
   field.RegisterParticles(&swarm);
   swarm.RegisterField(&field);

   cout << "Field and particle registered.\n";
#endif

	// Output initial data to a file
	FILE *outFile;
	if (argc==1)
	{
		outFile = fopen("debug.dat","a");
	}
	else
	{
		outFile = fopen(argv[1],"a");
	}

	for (int snaps=0; snaps < snapshots; snaps++)
	{
		// Relax grid using RK2 until convergence
		cout << "Starting relaxation routine for orbital separation " << simParams.particleCenter1.Val(1)-simParams.particleCenter0.Val(1) << ". \n";
		cout << "Current frame rotation rate is " << simParams.omega << "\n";

		simParams.omega = field.relax();

		cout << "Grid is relaxed, dumping initial data to output file.\n";

		for (i = 1; i<(int)simParams.max.Val(1)-1; i++)
		{
			for (j = 1; j<(int)simParams.max.Val(2)-1; j++)
			{
				for (k = 1; k<(int)simParams.max.Val(3)-1; k++)
				{
					fprintf(outFile,"%f %f %f %f %f %f",0.,
						field.Pos(1,i,j,k),
						field.Pos(2,i,j,k),
						field.Pos(3,i,j,k),
						field.Phi(1,i,j,k),
						field.Pi(1,i,j,k));
					fprintf(outFile,"\n");
				}
			}
		}

		fprintf(outFile,"%f %f %e %e",2.*simParams.particleCenter0.Val(1),simParams.omega,field.calcDipolePoyntingFlux(),field.calcQuadPoyntingFlux());
		fprintf(outFile,"\n");
		fflush(outFile);

		// Decrease the grid spacing, which has the effect of moving the particles inward
		gridspacing -= deltaspacing; 
		simParams.delta.Set(gridspacing/4.,0);
		simParams.delta.Set(gridspacing,1);
		simParams.delta.Set(gridspacing,2);
		simParams.delta.Set(gridspacing,3);

		// Calculate the particles' current separation
		R = nR*gridspacing;
		simParams.particleCenter0.Set(-R,1);
		simParams.particleCenter1.Set(R,1);
	
		// Set the radius of extraction
		simParams.re = nRe*gridspacing;

		// Transfer the new parameters to the field object
		field.SetParams(simParams);

	}

	// Close all files
	fclose(outFile);

	return(0);
}
