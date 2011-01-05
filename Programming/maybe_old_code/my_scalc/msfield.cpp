// MSField.cpp
// Implementation of Misner scalar field class
// Keith Watt 1999

#include "msfield.hpp"


MSField::MSField()
{
#ifdef PART
   my_particle_partner = NULL;
   matter_piece.Resize(2,4,4);
   field_piece. Resize(2,4,4);
#endif

   g.           Resize(2,4,4);
   dgdphi.      Resize(2,4,4);
   g_inv.       Resize(2,4,4);
   four_pos.    Resize(1,4);

   return;
}



void MSField::SetParams(Param simParamsIn)
{
	int q,i,j,k;

	simParams <= simParamsIn;
	omega = simParams.omega;

	//allocate an array of double pointers to doubles 
	phi = new double*** [3];
	pi  = new double*** [3];

	//now allocate an array of pointers to doubles for each member
	//of the array of double pointers to doubles
	//
	for(q = 0; q < 3; q++)
	{
		phi[q] = new double** [(int)simParams.max.Val(1)];
		pi[q]  = new double** [(int)simParams.max.Val(1)];
	}

	//allocate an array of doubles for each member of the arrays
	//of pointers to doubles
	for(q = 0; q < 3; q++)
	{
		for(i = 0; i < simParams.max.Val(1); i++)
		{
			phi[q][i] = new double* [(int)simParams.max.Val(2)];
			pi[q][i]  = new double* [(int)simParams.max.Val(2)];
		}
	}

	//allocate an array of doubles for each member of the arrays
	//of pointers to doubles
	for(q = 0; q < 3; q++)
	{
		for(i = 0; i < simParams.max.Val(1); i++)
		{
			for(j = 0; j < simParams.max.Val(2); j++)
			{
				phi[q][i][j] = new double [(int)simParams.max.Val(3)];
				pi[q][i][j]  = new double [(int)simParams.max.Val(3)];
			}
		}
	}
 
	//finally lets set the value of the array to be equal to the 1.0
	for (q = 0; q < 3; q++)
	{
		for(i = 0; i < simParams.max.Val(1); i++)
		{
			for( j = 0; j < simParams.max.Val(2); j++)
			{
				for(k = 0; k < simParams.max.Val(3); k++)
				{
					phi[q][i][j][k] = 0.0;
					pi[q][i][j][k]  = 0.0;
				}
			}
		}
	}

#ifdef skippy  //conditional compilation for debug purposes (randomize field)
  //set the seed
  srand(19);
  //guarentee that the random number lie between zero and -1
  double rand_norm = -1.0 / 32767.0;
  for (int i = 0; i < phi.p->product; i++)
  {
    phi.p->m[i] =  rand_norm * rand();
    pi.p->m[i]  =  rand_norm * rand();
  }
/*
  char *dump;

  FILE *my_dump_file;
  my_dump_file = fopen("init_phi.txt","w");

  phi.print("phi",&dump);
  fprintf(my_dump_file,"%s",dump);

  fclose(my_dump_file);
*/
  cout << "Debug checks:  phi[0,1,2,3] = "  << phi.Val(0,1,2,3)  << endl;
  cout << "Debug checks:  phi[2,10,6,5] = " << phi.Val(2,10,6,5) << endl;
  cout << "Debug checks:  phi[1,3,7,12] = " << phi.Val(1,3,7,12) << endl;
#endif

}



tensor MSField::Pos(int i, int j, int k)
{
  // Return a tensor containing the coordinate position of a grid-indexed
  // location
  // There is an inline version which returns only a single component
  // defined in MSField.hpp

  // Set the time value of the returned 4-vector to 0.  I (KW) always use
  // 4-vectors so that x is always index 1, y is always index 2, and
  //z is always index 3

  //temp.Set(0.,0);
  //commented out by CS - since temp (which is now rename four_pos)
  //is already initialized to be zero one need to nothing here

  // These are formulae I (KW) derived to calculate x, y, and z based upon the
  // grid size and spacing.
  // Note that the grid maximum includes the boundary cells
  //(unlike the Java-based ScalGrav34)
	four_pos.Set(((double)(i)-0.5*(simParams.max.Val(1)-1.0))*simParams.delta.Val(1),1);
	four_pos.Set(((double)(j)-0.5*(simParams.max.Val(2)-1.0))*simParams.delta.Val(2),2);
	four_pos.Set(((double)(k)-0.5*(simParams.max.Val(3)-1.0))*simParams.delta.Val(3),3);

	return four_pos;
}



// Find the closest grid location to a specified coordinate position
//	NOTE: This actually doesn't necessarily find the closest since it just rounds
//		all the numbers down.  It should probably check which of the bordering grid
//		locations is actually the closest.
int MSField::Index(int component, double x)
{
	// These are formulae I derived to caluclate x, y, and z based upon the grid size and spacing.
	// Note that the grid maximum includes the boundary cells (unlike the Java-based ScalGrav34)
	int ind = 0;

	// These formulae are essentially the inverse of the Pos formulae.  Again note that in this
	//	version the grid maximums include the boundary, unlike the Java-based ScalGrav34
	if (component == 1)
	{
		ind =  (int)(0.5*(simParams.max.Val(1)-1.) + x/simParams.delta.Val(1));
	}
	else if(component == 2)
	{
		ind =  (int)(0.5*(simParams.max.Val(2)-1.) + x/simParams.delta.Val(2));
	}
	else if(component == 3)
	{
		ind =  (int)(0.5*(simParams.max.Val(3)-1.) + x/simParams.delta.Val(3));
	}
	else
	{
		cout << "ERROR: Cannot find position of coordinate "<< component << "\n";
	}

	return ind;
}



// The entire definition of dPhi/dt
double MSField::dPhidt(double stepFrac, int q, int i, int j, int k)
{
	double x = Pos(1,i,j,k);
	double y = Pos(2,i,j,k);

	double dphidx = (phi[q][i+1][j][k]-phi[q][i-1][j][k])/(2.0*simParams.delta.Val(1));
	double dphidy = (phi[q][i][j+1][k]-phi[q][i][j-1][k])/(2.0*simParams.delta.Val(2));

	return stepFrac*simParams.delta.Val(0)*(pi[q][i][j][k]
			+ simParams.omega*(-y*dphidx + x*dphidy));
}



// The entire definition of dPi/dt
double MSField::dPidt(double stepFrac, double d, int q, int i, int j, int k)
{
	double source = 0.;

	if (simParams.dynCase == 0)
	{
		source = calcStaticSource(q,i,j,k);
	}
	else
	{
#ifdef PART
		matter_piece = my_particle_partner->matterVar(q,i,j,k);
		field_piece  = calcDgDphi(q,i,j,k);
		source = ( matter_piece.Contract(field_piece,1,0) ).Contract(0,1).Val(0);
#endif
	}


	source = -4.0*PI*G* source;
	double x = Pos(1,i,j,k);
	double y = Pos(2,i,j,k);

	double dphidx = (phi[q][i+1][j][k]-phi[q][i-1][j][k])/(2.0*simParams.delta.Val(1));
	double dphidy = (phi[q][i][j+1][k]-phi[q][i][j-1][k])/(2.0*simParams.delta.Val(2));
	double dpidx  = (pi[q][i+1][j][k]-pi[q][i-1][j][k])/(2.0*simParams.delta.Val(1));
	double dpidy  = (pi[q][i][j+1][k]-pi[q][i][j-1][k])/(2.0*simParams.delta.Val(2));
	double del2phi= (phi[q][i+1][j][k] - 2.0*phi[q][i][j][k] + phi[q][i-1][j][k])
						/pow(simParams.delta.Val(1),2.0)
					+(phi[q][i][j+1][k] - 2.0*phi[q][i][j][k] + phi[q][i][j-1][k])
						/pow(simParams.delta.Val(2),2.0)
					+(phi[q][i][j][k+1] - 2.0*phi[q][i][j][k] + phi[q][i][j][k-1])
						/pow(simParams.delta.Val(3),2.0);


	return stepFrac*simParams.delta.Val(0)
				*(del2phi + simParams.omega*(-y*dpidx + x*dpidy)
				-2.0*d*pi[q][i][j][k] - 2.0*d*simParams.omega*(-y*dphidx+x*dphidy)
				+ source);
}

// Calculate the total source term for the static case
double MSField::calcStaticSource(int q, int x, int y, int z)
{

	// q=timeslice (old, current, or new)
	// x,y,z=grid location

	double totalSource = 0.;
	// Add up the contribution from each particle
	for (int a=0; a<2; a++)
	{
		// Calculate the squared distance from the particle's center
		//	to the grid location specified as a fraction of the
		//	smoothing length squared
		double pcx,pcy,pcz,sl,mass;
		int    pci,pcj,pck;

		if (a==0)
		{
			pcx = simParams.particleCenter0.Val(1);
			pcy = simParams.particleCenter0.Val(2);
			pcz = simParams.particleCenter0.Val(3);
			sl  = simParams.smoothLength0.Val(1);
			mass = simParams.mass0;
		}
		else if (a==1)
		{
			pcx = simParams.particleCenter1.Val(1);
			pcy = simParams.particleCenter1.Val(2);
			pcz = simParams.particleCenter1.Val(3);
			sl  = simParams.smoothLength1.Val(1);
			mass = simParams.mass1;
		}

		pci = Index(1,pcx);
		pcj = Index(2,pcy);
		pck = Index(3,pcz);

		double r2 = (pow(pcx-Pos(1,x,y,z),2.)
				+ pow(pcy-Pos(2,x,y,z),2.)
				+ pow(pcz-Pos(3,x,y,z),2.))
			/pow(sl,2.0);

		// w is the weighting kernel
		double w = 0.0;

		// norm is the normalization of the weighting kernel
		double norm = 1.0;
		if (r2 < 1.0)
		{
			// Calculate normalization constant	by "smoothing" a value of 1
			norm = 0.0;
			double rnorm2;

			// It's not necessary to scan the entire grid, since
			//	anything outside the smoothing length is zero.
			// NOTE: It would probably be a good idea to add some
			//	"Tennessee windage" to this since the Index routine
			//	may not reliably giving the closest grid location.
			int minX = 0; //Index(1,pcx-sl*1.5);
			int maxX = (int)simParams.max.Val(1)-1; //Index(1,pcx+sl*1.5);
			int minY = 0; //Index(2,pcy-sl*1.5);
			int maxY = (int)simParams.max.Val(2)-1; //Index(2,pcy+sl*1.5);
			int minZ = 0; //Index(3,pcz-sl*1.5);
			int maxZ = (int)simParams.max.Val(3)-1; //Index(3,pcz+sl*1.5);

			// Scan within the mins and maxs defined above
			for (int i=minX; i<=maxX; i++)
			{
				for (int j=minY; j<=maxY; j++)
				{
					for (int k=minZ; k<=maxZ; k++)
					{
						// Calculate the squared distance between the particle
						//	center and the grid location being scanned as a
						//	fraction of the smoothing length squared

						rnorm2 = (pow(pcx-Pos(1,i,j,k),2.)
								+ pow(pcy-Pos(2,i,j,k),2.)
								+ pow(pcz-Pos(3,i,j,k),2.))
							/pow(sl,2.0);

						if (rnorm2 < 1.0)
						{
							// The smoothing kernel Charlie and I worked out
							//	Note that we are still "smoothing" a value of 1 here
							norm  = norm + pow(1.0-rnorm2,simParams.exponent);
						}
					}
				}
			}
			// Calculate the actual weighting of the requested grid location
			w = pow(1.0-r2,simParams.exponent);
		}

		// Check to see if omega is too large
		if (simParams.omega >=fabs(1./(pcx*exp(-2.*phi[1][pci][pcj][pck]))))
		{
			simParams.omega =fabs(1./(pcx*exp(-2.*phi[1][pci][pcj][pck])) * 0.99);
			cout << "WARNING: Calculated omega too large for static source, adjusting to maximum omega: " << simParams.omega << "\n";
		}

		if (w!=0.)
		{
		// Calculate the contribution of this particle to the total static source term
		double k = simParams.delta.Val(1);
		totalSource += w/(norm*k*k*k)*mass
	        * exp(phi[q][pci][pcj][pck]) * (1.+ simParams.omega*simParams.omega*pcx*pcx*exp(-4.*phi[q][pci][pcj][pck]))
				/ sqrt(1.-simParams.omega*simParams.omega*pcx*pcx*exp(-4.*phi[q][pci][pcj][pck]));
		}
	}
	
	return totalSource;
}

// Calculate the bare mass of the particles on the grid
double MSField::calcBareMass()
{

	double bareMass = 0.;

	// Add up the contribution from each particle
	for (int a=0; a<2; a++)
	{

		// Calculate the squared distance from the particle's center
		//	to the grid location specified as a fraction of the
		//	smoothing length squared
		double pcx,pcy,pcz,sl,mass;

		if (a==0)
		{
			pcx = simParams.particleCenter0.Val(1);
			pcy = simParams.particleCenter0.Val(2);
			pcz = simParams.particleCenter0.Val(3);
			sl  = simParams.smoothLength0.Val(1);
			mass = simParams.mass0;
		}
		else if (a==1)
		{
			pcx = simParams.particleCenter1.Val(1);
			pcy = simParams.particleCenter1.Val(2);
			pcz = simParams.particleCenter1.Val(3);
			sl  = simParams.smoothLength1.Val(1);
			mass = simParams.mass1;
		}

		// w is the weighting kernel
		double w = 0.;

		// norm is the normalization of the weighting kernel
		double norm = 0.;

		// Calculate normalization constant	by "smoothing" a value of 1
		double rnorm2;

		// It's not necessary to scan the entire grid, since
		//	anything outside the smoothing length is zero.
		// NOTE: It would probably be a good idea to add some
		//	"Tennessee windage" to this since the Index routine
		//	may not reliably giving the closest grid location.
		int minX = 0; //Index(1,pcx-sl*1.5);
		int maxX = (int)simParams.max.Val(1)-1; //Index(1,pcx+sl*1.5);
		int minY = 0; //Index(2,pcy-sl*1.5);
		int maxY = (int)simParams.max.Val(2)-1; //Index(2,pcy+sl*1.5);
		int minZ = 0; //Index(3,pcz-sl*1.5);
		int maxZ = (int)simParams.max.Val(3)-1; //Index(3,pcz+sl*1.5);

		// Scan within the mins and maxs defined above
		for (int i=minX; i<=maxX; i++)
		{
			for (int j=minY; j<=maxY; j++)
			{
				for (int k=minZ; k<=maxZ; k++)
				{
					// Calculate the squared distance between the particle
					//	center and the grid location being scanned as a
					//	fraction of the smoothing length squared

					rnorm2 = (pow(pcx-Pos(1,i,j,k),2.)
							+ pow(pcy-Pos(2,i,j,k),2.)
							+ pow(pcz-Pos(3,i,j,k),2.))
						/pow(sl,2.0);

					if (rnorm2 < 1.0)
					{
						// The smoothing kernel Charlie and I worked out
						//	Note that we are still "smoothing" a value of 1 here
						norm  += pow(1.0-rnorm2,simParams.exponent);
					}
				}
			}
		}

		// Scan within the mins and maxs defined above
		for (int x=minX; x<=maxX; x++)
		{
			for (int y=minY; y<=maxY; y++)
			{
				for (int z=minZ; z<=maxZ; z++)
				{

					double r2 = (pow(pcx-Pos(1,x,y,z),2.)
							+ pow(pcy-Pos(2,x,y,z),2.)
							+ pow(pcz-Pos(3,x,y,z),2.))
						/pow(sl,2.0);

					if (r2 < 1.0)
					{
						// Calculate the actual weighting of the requested grid location
						w = pow(1.0-r2,simParams.exponent);
					}
					else
					{
						w = 0.;
					}

					bareMass += mass*w/norm;
				}
			}
		}
	}
	
	return bareMass;
}

double MSField::calcOmega(double currOmega)
{

	int pcx = Index(1,simParams.particleCenter0.Val(1));
	int pcy = Index(2,simParams.particleCenter0.Val(2));
	int pcz = Index(3,simParams.particleCenter0.Val(3));
	double r = fabs(simParams.particleCenter0.Val(1));


	double dphidx = (phi[1][pcx+1][pcy][pcz]-phi[1][pcx-1][pcy][pcz])/(2.*simParams.delta.Val(1));
	double output = sqrt(fabs(dphidx/(r*exp(-4.*phi[1][pcx][pcy][pcz])*(1.-r*dphidx))));


	if (output == 0.)
	{
		output = fabs(1./(r*exp(-2.*phi[1][pcx][pcy][pcz])) * 0.5);
		cout << "WARNING: Calculated omega equals zero, adjusting to maximum omega. \n";
	}

//	cout << "Calculation complete. \n";

	return output;
}

double MSField::calcMotorForce(int i)
{
	int pcx = Index(1,simParams.particleCenter0.Val(1));
	int pcy = Index(2,simParams.particleCenter0.Val(2));
	int pcz = Index(3,simParams.particleCenter0.Val(3));
	double r = fabs(simParams.particleCenter0.Val(1));
	double output = 0.;

	double dphidi = 0.;

	if (i==1)
	{
		dphidi = (phi[1][pcx+1][pcy][pcz]-phi[1][pcx-1][pcy][pcz])/(2.*simParams.delta.Val(1));
	}
	else if (i==2)
	{
		dphidi = (phi[1][pcx][pcy+1][pcz]-phi[1][pcx][pcy-1][pcz])/(2.*simParams.delta.Val(2));
	}
	else if (i==3)
	{
		dphidi = (phi[1][pcx][pcy][pcz+1]-phi[1][pcx][pcy][pcz-1])/(2.*simParams.delta.Val(3));
	}

	output = exp(phi[1][pcx][pcy][pcz])*dphidi*(1.+simParams.omega*simParams.omega*r*r*exp(-4.*phi[1][pcx][pcy][pcz]))
				/ sqrt(1.-simParams.omega*simParams.omega*r*r*exp(-4.*phi[1][pcx][pcy][pcz]));

	if (i==1)
	  {
	    output += -r*simParams.omega*simParams.omega
			*exp(-4.*phi[1][pcx][pcy][pcz])/sqrt(1.-simParams.omega*simParams.omega*r*r*exp(-4.*phi[1][pcx][pcy][pcz]));
	  }

	return output;
}

double MSField::calcDgDphi(int q, int i, int j, int k,  int a, int b)
{
	double dgdphi;
	double x = Pos(1,i,j,k);
	double y = Pos(2,i,j,k);
	double r2 = x*x + y*y;

	if (a==0 && b==0)
	{
		dgdphi = 2.*exp(-2.*phi[q][i][j][k]);
	}
	else if ((a==0 && b==1) || (a==1 && b==0))
	{
		dgdphi = 2*simParams.omega*y*exp(-2.*phi[q][i][j][k]);
	}
	else if ((a==0 && b==2) || (a==2 && b==0))
	{
		dgdphi = 2.*simParams.omega*x*exp(-2.*phi[q][i][j][k]);
	}
	else if ((a==0 && b==3) || (a==3 && b==0))
	{
		dgdphi = 0.;
	}
	else if (a==1 && b==1)
	{
		dgdphi = 2.*exp(2.*phi[q][i][j][k])
					+simParams.omega*simParams.omega*y*y*exp(-2.*phi[q][i][j][k]);
	}
	else if ((a==1 && b==2) || (a==2 && b==1))
	{
		dgdphi = -2.*x*y*simParams.omega*simParams.omega*exp(-2.*phi[q][i][j][k]);
	}
	else if ((a==1 && b==3) || (a==1 && b==3))
	{
		dgdphi = 0.;
	}
	else if (a==2 && b==2)
	{
		dgdphi = 2.*exp(2.*phi[q][i][j][k])
					+simParams.omega*simParams.omega*x*x*exp(-2.*phi[q][i][j][k]);
	}
	else if ((a==2 && b==3) || (a==3 && b==2))
	{
		dgdphi = 0.;
	}
	else if ((a==3 && b==3) || (a==3 && b==3))
	{
		dgdphi = 2.*exp(2.*phi[q][i][j][k]);
	}

	return dgdphi;
}


//*****************************************************************************
//Name:   calcDgDphi
//
//Purpose: returns the result of the partial of the inverse metric with
//         respect to the scalar field phi
//
//Takes:  time slice and position labels as ints
//
//Written by:	CS 2/5/00
//******************************************************************************
tensor MSField::calcDgDphi(int q, int i, int j, int k)
{
   double phi_value;
   double exp_m2_phi;
   double exp_p2_phi;

   double x = Pos(1,i,j,k);
   double y = Pos(2,i,j,k);
   double r2 = x*x + y*y;

   phi_value  = phi[q][i][j][k];
   exp_m2_phi = exp(-2.0*phi_value);
   exp_p2_phi = exp(+2.0*phi_value);

   dgdphi.p->m[0]  = 2.0*exp_m2_phi;                    //0 0 component
   dgdphi.p->m[1]  = 2.0*omega*y*exp_m2_phi;            //0 1 component
   dgdphi.p->m[2]  = -2.0*omega*x*exp_m2_phi;           //0 2 component
   dgdphi.p->m[3]  = 0.0;                               //0 3 component
   dgdphi.p->m[4]  = 2.0*omega*y*exp_m2_phi;            //1 0 component
   dgdphi.p->m[5]  = 2.0*exp_p2_phi +
                        2.0*omega*omega*y*y*exp_m2_phi; //1 1 component
   dgdphi.p->m[6]  = -2.0*x*y*omega*omega*exp_m2_phi;   //1 2 component
   dgdphi.p->m[7]  = 0.0;                               //1 3 component
   dgdphi.p->m[8]  = -2.0*omega*x*exp_m2_phi;           //2 0 component
   dgdphi.p->m[9]  = -2.0*x*y*omega*omega*exp_m2_phi;   //2 1 component
   dgdphi.p->m[10] = 2.0*exp_p2_phi +
                        2.0*omega*omega*x*x*exp_m2_phi; //2 2 component
   dgdphi.p->m[11] = 0.0;                               //2 3 component
   dgdphi.p->m[12] = 0.0;                               //3 0 component
   dgdphi.p->m[13] = 0.0;                               //3 1 component
   dgdphi.p->m[14] = 0.0;                               //3 2 component
   dgdphi.p->m[15] = 2.0*exp_p2_phi;                    //3 3 component

   return dgdphi;
}


// Return the contravariant metric tensor on the grid
// q is the time slice (0="old", 1="current", 2="new")
// i,j,k is the grid coordinate
double MSField::metricG(int q, int i, int j, int k, int a, int b)
{
	double g;
	double x = Pos(1,i,j,k);
	double y = Pos(2,i,j,k);
	double r2 = x*x + y*y;

	if (a==0 && b==0)
	{
		g = -exp(2.*phi[q][i][j][k]) + r2*simParams.omega*simParams.omega*exp(-2.*phi[q][i][j][k]);
	}
	else if ((a==0 && b==1) || (a==1 && b==0))
	{
		g = -y*simParams.omega*exp(-2.*phi[q][i][j][k]);
	}
	else if ((a==0 && b==2) || (a==2 && b==0))
	{
		g = x*simParams.omega*exp(-2.*phi[q][i][j][k]);
	}
	else if ((a==0 && b==3) || (a==3 && b==0))
	{
		g = 0.;
	}
	else if (a==1 && b==1)
	{
		g = exp(-2.*phi[q][i][j][k]);
	}
	else if ((a==1 && b==2) || (a==2 && b==1))
	{
		g = 0.;
	}
	else if ((a==1 && b==3) || (a==1 && b==3))
	{
		g = 0.;
	}
	else if (a==2 && b==2)
	{
		g = exp(-2.*phi[q][i][j][k]);
	}
	else if ((a==2 && b==3) || (a==3 && b==2))
	{
		g = 0.;
	}
	else if ((a==3 && b==3) || (a==3 && b==3))
	{	
		g = exp(-2.*phi[q][i][j][k]);
	}

	return g;
}

// Return the contravariant metric tensor on the grid
// q is the time slice (0="old", 1="current", 2="new")
// i,j,k is the grid coordinate
//updated by CS 2/9/00
tensor MSField::GetInvG(int q, int i, int j, int k)
{
   double phi_value;
   double exp_m2_phi;
   double exp_p2_phi;

	double x = Pos(1,i,j,k);
	double y = Pos(2,i,j,k);
	double r2 = x*x + y*y;

   phi_value  = phi[q][i][j][k];
   exp_m2_phi = exp(-2.0*phi_value);
   exp_p2_phi = exp(+2.0*phi_value);

   g_inv.p->m[0]  = -exp_m2_phi;                       //0 0 component
   g_inv.p->m[1]  = -omega*y*exp_m2_phi;               //0 1 component
   g_inv.p->m[2]  = omega*x*exp_m2_phi;                //0 2 component
   g_inv.p->m[3]  = 0.0;                               //0 3 component
   g_inv.p->m[4]  = -omega*y*exp_m2_phi;               //1 0 component
   g_inv.p->m[5]  = exp_p2_phi -
                        omega*omega*y*y*exp_m2_phi;    //1 1 component
   g_inv.p->m[6]  = x*y*omega*omega*exp_m2_phi;        //1 2 component
   g_inv.p->m[7]  = 0.0;                               //1 3 component
   g_inv.p->m[8]  = omega*x*exp_m2_phi;                //2 0 component
   g_inv.p->m[9]  = x*y*omega*omega*exp_m2_phi;        //2 1 component
   g_inv.p->m[10] = exp_p2_phi -
                        omega*omega*x*x*exp_m2_phi;    //2 2 component
   g_inv.p->m[11] = 0.0;                               //2 3 component
   g_inv.p->m[12] = 0.0;                               //3 0 component
   g_inv.p->m[13] = 0.0;                               //3 1 component
   g_inv.p->m[14] = 0.0;                               //3 2 component
   g_inv.p->m[15] = exp_p2_phi;                        //3 3 component

	return g_inv;
/*
	tensor g(2,4,4);
	double x = Pos(1,i,j,k);
	double y = Pos(2,i,j,k);

	g.Set(-exp(-2.*phi.Val(q,i,j,k)),0,0);
	g.Set(-simParams.omega*y*exp(-2.*phi.Val(q,i,j,k)),0,1);
	g.Set(simParams.omega*x*exp(-2.*phi.Val(q,i,j,k)),0,2);
	g.Set(0.,0,3);

	g.Set(-simParams.omega*y*exp(-2.*phi.Val(q,i,j,k)),1,0);
	g.Set(exp(2.*phi.Val(q,i,j,k)) - y*y*simParams.omega*simParams.omega*exp(-2.*phi.Val(q,i,j,k)),1,1);
	g.Set(x*y*simParams.omega*simParams.omega*exp(-2.*phi.Val(q,i,j,k)),1,2);
	g.Set(0.,1,3);

	g.Set(simParams.omega*x*exp(-2.*phi.Val(q,i,j,k)),2,0);
	g.Set(x*y*simParams.omega*simParams.omega*exp(-2.*phi.Val(q,i,j,k)),2,1);
	g.Set(exp(2.*phi.Val(q,i,j,k)) - x*x*simParams.omega*simParams.omega*exp(-2.*phi.Val(q,i,j,k)),2,2);
	g.Set(0.,2,3);

	g.Set(0.,3,0);
	g.Set(0.,3,1);
	g.Set(0.,3,2);
	g.Set(exp(2.*phi.Val(q,i,j,k)),3,3);

	return g;
*/
}


// Calculate the energy of the grid based upon the Lagrangian of the
//	scalar gravity field equation
double MSField::calcEnergy(int q)
{
	/* PI and G are defined in constant.hpp
		NOTE: Presently this routine ignores the contribution due to the matter
		Integrate the Lagrangian over the entire grid.  Note that I 
		don't presently include the boundaries in the intergration.
	*/
	tensor t(2,4,4);

	double energy = 0.0;
	for (int i = 1; i < (int)simParams.max.Val(1)-1; i++)
	{
		 for (int j = 1; j < (int)simParams.max.Val(2)-1; j++)
		 {
			for (int k = 1; k < (int)simParams.max.Val(3)-1; k++)
			{
				 energy += 0.5*pow(pi[q][i][j][k],2.0)
					+ 0.5*pow(
					 (phi[q][i+1][j][k] - phi[q][i-1][j][k])/(2.0*simParams.delta.Val(1))
					+(phi[q][i][j+1][k] - phi[q][i][j-1][k])/(2.0*simParams.delta.Val(2))
					+(phi[q][i][j][k+1] - phi[q][i][j][k-1])/(2.0*simParams.delta.Val(3))
					,2.0);
			}
		 }
	}

	return energy
		*simParams.delta.Val(1)*simParams.delta.Val(2)*simParams.delta.Val(3);
}

// Fill the grid with zero values of phi and pi
void MSField::flatSpace()
{
	for (int q=0; q<3; q++)
	{
		for (int i=0; i<(int)simParams.max.Val(1); i++)
		{
			for (int j=0; j<(int)simParams.max.Val(2); j++)
			{
				for (int k=0; k<(int)simParams.max.Val(3); k++)
				{
					phi[q][i][j][k] = 0.;
					pi[q][i][j][k]  = 0.;
				}
			}
		}
	}

	return;
}


// Relax the grid by integrating with diffusion.  When the energy of the grid is unchanging
//	by more than the convergence goal, we'll consider the grid to have relaxed.
double MSField::relax()
{

	// Initialize the harmonic decomposition if necessary
	if (simParams.boundCond == HARM_DECOMP)
	{
		cout << "Initializing harmonic decomposition.\n";
		initMisnerHarmonics();
	}

	// itNum is the iteration counter, used to give an idea of how many steps are required
	//	for convergence (typically ~200)
	int itNum = 0;
	// itCount is used to track the number of iterations since the grid energy was lasted output
	//	to the screen.
	int itCount = 0;
	// remainder stores essentially the time derivative of the grid energy.  Ideally, this should be zero.
	double remainder;

	// currEnergy is the current grid energy
	double currEnergy = 0.;
	// oldEnergy is the energy of the grid from the previous iteration.  Used in calculating the time
	//	derivative for remainder
	double oldEnergy;

	// Open an output file for monitor data
	FILE *monitorFile;
	monitorFile = fopen("monitor.dat","a");

	double bareMass = calcBareMass();

	cout << "Initial Bare Mass = " << bareMass << "\n";
	fprintf(monitorFile,"%f",bareMass);

	do
	{
		// Take an RK2 step to evolve the grid
		// NOTE: It's not required that the relaxation timestep be the same as the 
		//	system timestep, but I don't see any reason to make it otherwise
		rk2(1.0, simParams.diffCoeff);

		// Store the current energy in oldEnergy and calculate the new energy
		oldEnergy = currEnergy;
		currEnergy = calcQuadPoyntingFlux(); //calcEnergy(1);

		// Calculate the time derivative of the energy.
		// The cases below are to avoid having divide by zero errors.
		if (currEnergy != 0.)
		{
			remainder = (currEnergy - oldEnergy)/(currEnergy*simParams.delta.Val(0));
		}
		else if (oldEnergy != 0.)
		{
			remainder = 1.0/simParams.delta.Val(0);
		}
		else
		{
			remainder = 100.0;
		}

		itNum++;
		itCount++;

		// NOTE: This is currently hard-coded to print every iteration, but an equivalent to
		//	printSkip (or even printSkip itself) could be used here.
		if (itCount==1)
		{
			itCount = 0;

			int pcx = Index(1,simParams.particleCenter0.Val(1));
			int pcy = Index(2,simParams.particleCenter0.Val(2));
			int pcz = Index(3,simParams.particleCenter0.Val(3));
			double dphidx = (phi[1][pcx+1][pcy][pcz]-phi[1][pcx-1][pcy][pcz])/(2.*simParams.delta.Val(1));
			double ph = phi[1][pcx][pcy][pcz];

			// Print the derivative and the current grid energy
			cout << simParams.particleCenter1.Val(1) << " " << itNum << " " << remainder << " " << currEnergy << " " << simParams.omega << " " << dphidx << " " << ph << "\n";
			fprintf(monitorFile,"%f %i %e %e %f %e %f",simParams.particleCenter1.Val(1),itNum,remainder,currEnergy,simParams.omega,dphidx,ph);
			fprintf(monitorFile,"\n");
			fflush(monitorFile);

		}
		simParams.omega = calcOmega(simParams.omega);
  }
  // Loop until we meet the convergence goal
  while(fabs(remainder)>simParams.convGoal);
	
  for (int a = 0; a<numHarm(); a++)
  {
		cout << "b" << a << "=" << b[a] << "\n";
		fprintf(monitorFile, "%e", b[a]);
		fprintf(monitorFile,"\n");
  }

	cout << "Dipole Poynting Flux = " << calcDipolePoyntingFlux() << "\n";
	fprintf(monitorFile,"%e %e",calcDipolePoyntingFlux(), calcQuadPoyntingFlux());
	fprintf(monitorFile,"\n");

	cout << "Quadrupole Poynting Flux = " << calcQuadPoyntingFlux() << "\n";

	cout << "Final Frame Rotation Rate = " << simParams.omega << "\n";
	fprintf(monitorFile,"%f",simParams.omega);
	fprintf(monitorFile,"\n");

	cout << "Final Motor x-Force = " << calcMotorForce(1) << "\n";
	cout << "Final Motor y-Force = " << calcMotorForce(2) << "\n";
	cout << "Final Motor z-Force = " << calcMotorForce(3) << "\n";
	fprintf(monitorFile,"%e %e %e",calcMotorForce(1),calcMotorForce(2),calcMotorForce(3));
	fprintf(monitorFile,"\n");

	bareMass = calcBareMass();

	cout << "Final Bare Mass = " << bareMass << "\n";
	fprintf(monitorFile,"%f",bareMass);

	fclose(monitorFile);


  return simParams.omega;
}


// Take an RK2 step; stepFrac is the fraction of the system timestep,
//	d is the diffusion coefficient used by the relaxation routine
void MSField::rk2(double stepFrac, double d)
{
	int i,j,k;

	// Scan the entire grid except for the boundaries
	for (i=1; i<(int)simParams.max.Val(1)-1; i++)
	{
		for (j=1; j<(int)simParams.max.Val(2)-1; j++)
		{
			for (k=1; k<(int)simParams.max.Val(3)-1; k++)
			{
				if (rDist(i,j,k) <= simParams.re + 1.5 * simParams.delta.Val(1))
				{
					// Find the midpoint (1/2*stepFrac of the system timestep) value of phi using the rotating coordinate equations
					phi[0][i][j][k] = phi[1][i][j][k] + dPhidt(0.5*stepFrac, 1,i,j,k);
					
					// Find the midpoint (1/2*stepFrac of the system timestep) value of pi using the rotating coordinate equations
					pi[0][i][j][k]  = pi[1][i][j][k]  + dPidt(0.5*stepFrac, d, 1,i,j,k);
				}
			}
		}
	}


	// Calculate boundary conditions; the 0 means to use the "old" (which is where RK2 stores its midpoint
	//	values) timeslice

	calcBoundary(0);

	// Second half of RK step
	// Scan the entire grid except for the boundaries
	for (i=1; i<(int)simParams.max.Val(1)-1; i++)
	{
		for (j=1; j<(int)simParams.max.Val(2)-1; j++)
		{
			for (k=1; k<(int)simParams.max.Val(3)-1; k++)
			{
				if (rDist(i,j,k) <= simParams.re + 1.5 * simParams.delta.Val(1))
				{
					// Find the new (stepFrac * the system timestep) value of phi using the rotating coordinate equations
					phi[2][i][j][k] = phi[1][i][j][k] + dPhidt(stepFrac, 0,i,j,k);

					// Find the midpoint (new stepFrac * the system timestep) value of pi using the rotating coordinate equations
					pi[2][i][j][k]  = pi[1][i][j][k] + dPidt(stepFrac, d, 0,i,j,k);
				}
			}
		}
	}

	// Calculate boundary conditions; 2 means to use the "new" timeslice
	calcBoundary(2);
	
	// Update the field variables, moving the current into old
	//		and new into current
	for (i=0; i<(int)simParams.max.Val(1); i++)
	{
		for (j=0; j<(int)simParams.max.Val(2); j++)
		{
			for (k=0; k<(int)simParams.max.Val(3); k++)
			{
				phi[0][i][j][k] = phi[1][i][j][k];
				pi[0][i][j][k]  = pi[1][i][j][k];
				phi[1][i][j][k] = phi[2][i][j][k];
				pi[1][i][j][k]  = pi[2][i][j][k];
			}
		}
	}

	return;
}

// Calculate the boundary for the timeslice q
//	0 = "old" (midpoint for RK2)
//	1 = "current"
//	2 = "new"
void MSField::calcBoundary(int q)
{
	int i,j,k;

	// Set the boundary values to zero
	// NOTE: This is hard-coded to be zero, but actually could be any fixed value
	if (simParams.boundCond==ZERO_FIXED)
	{
		// The +/- yz-face of the grid "cube"
		for (j = 0; j <= (int)simParams.max.Val(2)-1; j++)
		{
			for (k = 0; k <= (int)simParams.max.Val(3)-1; k++)
			{
				// Negative side
 				phi[q][0][j][k] = 0.;
				pi[q][0][j][k]  = 0.;
				// Positive side
				phi[q][(int)simParams.max.Val(1)-1][j][k] = 0.;
				pi[q][(int)simParams.max.Val(1)-1][j][k]  = 0.;
			}
		}

		// The +/- xz-face of the grid "cube"
		for (i = 0; i <= (int)simParams.max.Val(1)-1; i++)
		{
			for (k = 0; k <= (int)simParams.max.Val(3)-1; k++)
			{
				// Negative side
				phi[q][i][0][k] = 0.;
				pi[q][i][0][k]  = 0.;
				// Positive side
				phi[q][i][(int)simParams.max.Val(2)-1][k] = 0.;
				pi[q][i][(int)simParams.max.Val(2)-1][k]  = 0.;
			}
		}

		// The +/- xy-face of the grid "cube"
		for (i = 0; i <= (int)simParams.max.Val(1)-1; i++)
		{
			for (j = 0; j <= (int)simParams.max.Val(2)-1; j++)
			{
				// Negative side
				phi[q][i][j][0] = 0.;
				pi[q][i][j][0]  = 0.;
				// Positive side
				phi[q][i][j][(int)simParams.max.Val(3)-1] = 0.;
				pi[q][i][j][(int)simParams.max.Val(3)-1]  = 0.;
			}
		}
	}

	// Fill in the grid from the extraction radius outward with values
	//	calculated from the analytic solution to the field equation.
	// NOTE: This analytic solution is only valid where T=0, so
	//	the extraction radius must be outside of any matter, including 
	//	smoothing
	else if (simParams.boundCond == HARM_DECOMP)
	{
		int i,j,k;

		// Calculate new harmonic decomposition based on the current state of the grid.
		// This is necessary because the harmonic coefficients will change each timestep
		//	due to the non-linear portion of the grid
		calcMisnerHarmonics();

		// It's not necessary to scan the entire grid, since 
		//	we only fill outside the extraction radius, but I 
		//	haven't taken the time to work out the formulae
		//	to only scan the exterior regions.  
		for (i=0; i<=(int)simParams.max.Val(1)-1; i++)
		{
			for (j=0; j<=(int)simParams.max.Val(2)-1; j++)
			{
				for (k=0; k<=(int)simParams.max.Val(3)-1; k++)
				{
					if (rDist(i,j,k) > simParams.re)
					{
						phi[q][i][j][k] = restoreFromMisnerHarmonics(i,j,k);
					}
				}
			}
		}
	}

	// The specified boundary conditions haven't yet been implemented
	else
	{
		cout << "ERROR: Invalid boundary conditions\n";	// ERROR
	}
}

// Calculate the total number of harmonics present based upon the 
//	maximum value of l and the number of "Misner harmonic" (the
//	Legendre polynomial-based terms used for the grid-shell extraction)
//	terms to use
int MSField::numHarm()
{
	int numH = 0;
	for (int n=0; n<=simParams.maxL; n++)
	{
		// There are 2l+1 harmonics for each value of l
		numH += 2*n+1;
	}

	// Multiply the number of "regular" spherical harmonics
	//	by the number of "Misner" harmonics, e.g., if 3 terms in
	//	the Misner expansion are kept, there will be three for -each-
	//	spherical harmonic
	return numH*(simParams.maxN+1);
}

// We don't want to keep a ragged array of harmonic coefficients (since we are going to have
//	to invert the normalization tensor), so return a tensor which contains the harmonic indices
//	which correspond to the harmonic array index specified.
// NOTE: This routine significantly slows down the code.  A routine which returns only a specific 
//	harmonic index is implemented below.
tensor MSField::getK(int k)
{
	int n=-1;
	int m=0;
	int l=0;
		
	// Count up through the harmonic indices until we've counted a number of 
	//	harmonics equal to the requested array index
	for (int q=0; q<=k; q++)
	{
		n=n+1;
		if (n>simParams.maxN)
		{
			n=0;
			m=m+1;
			if (m>l)
			{
				l=l+1;
				m=-l;
			}
		}
	}
	
	// The harmonic indices are assigned to the tensor compoonents as follows:
	//	0 = l
	//	1 = m
	//	2 = n (Misner harmonic)
	tensor kVect(1,3);
	kVect.Set((double)l,0);
	kVect.Set((double)m,1);
	kVect.Set((double)n,2);
	
	return kVect;
}	

// We don't want to keep a ragged array of harmonic coefficients (since we are going to have
//	to invert the normalization tensor), so return the harmonic index which corresponds to the 
//	harmonic array index specified and the component requested.
int MSField::getK(int comp, int k)
{
	int n=-1;
	int m=0;
	int l=0;
		
	// Count up through the harmonic indices until we've counted a number of 
	//	harmonics equal to the requested array index
	for (int q=0; q<=k; q++)
	{
		n=n+1;
		if (n>simParams.maxN)
		{
			n=0;
			m=m+1;
			if (m>l)
			{
				l=l+1;
				m=-l;
			}
		}
	}
	
	if (comp==0)		// User wants the l-index
	{
		return l;
	}
	else if (comp==1)	// User wants the m-index
	{
		return m;
	}
	else if (comp==2)	// User wants the n-index
	{
		return n;
	}
	else				// Invalid index requested
	{
		cout << "ERROR: Attempt to calculate invalid harmonic index. \n";
		return 0;
	}
}	


// We don't want to keep a ragged array of harmonic coefficients (since we are going to have
//	to invert the normalization tensor), so return the harmonic array index which corresponds to the 
//	harmonic indices specified.
int MSField::getK(int l, int m, int n)
{
	int n1=-1;
	int m1=0;
	int l1=0;
	int k=-1;
		
	// Count up through the harmonic indices until we've reached the specified index set.
	//	Keep track of how many indices we had to count through - that will be the harmonic 
	//	array index
	do
	{
		k=k+1;
		n1=n1+1;
		if (n1>simParams.maxN)
		{
			n1=0;
			m1=m1+1;
			if (m1>l1)
			{
				l1=l1+1;
				m1=-l1;
			}
		}
	}
	while (l1!=l || m1!=m || n1!=n);
		
	return k;
}	

// Because we can't count on the numerical grid-shell extraction procedure yielding an 
//	orthogonal basis, calculate the harmonic decomposition normalization tensor.  This
//	tensor will be inverted and used when taking the harmonic decomposition at the 
//	extraction radius.	

tensor MSField::calcNormTensor()
{
	// Create a Util object, which has all the various special function (spherical harmonics,
	//	Legendre polynomials, etc.) we will need.
	Util util;

	// Resize the normaliztion tensor to the proper number of harmonic array indices
	g.Resize(2,numHarm(),numHarm());

	// It's not necessary to scan the grid since the Simpson's Rule weighting falls off
	//	inside and outside of the shell.  I've added a little "Tennessee windage" since the 
	//	Index routine doesn't necessarily give the actual closest grid location to a 
	//	coordinate position.
	int xMin = 0; //Index(1,-(simParams.re+9.*simParams.delta.Val(1)/4.));
	int xMax = (int)simParams.max.Val(1)-1; //Index(1,simParams.re+9.*simParams.delta.Val(1)/4.);
	int yMin = 0; //Index(2,-(simParams.re+9.*simParams.delta.Val(1)/4.));
	int yMax = (int)simParams.max.Val(2)-1; //Index(2,simParams.re+9.*simParams.delta.Val(1)/4.);
	int zMin = 0; //Index(3,-(simParams.re+9.*simParams.delta.Val(1)/4.));
	int zMax = (int)simParams.max.Val(3)-1; //Index(3,simParams.re+9.*simParams.delta.Val(1)/4.);
		
	// Scan within the shell
	for (int i = xMin; i<=xMax; i++)
	{
		for (int j = yMin; j<=yMax; j++)
		{
			for (int k = zMin; k<=zMax; k++)
			{
				double w = 0.0;

				// Calculate the coordinates of the grid location in spherical coordinates
				double r = rDist(i,j,k);
				if (r>=simParams.re-5.*simParams.delta.Val(1)/4. && r<=simParams.re+5.*simParams.delta.Val(1)/4.)
				{
					double theta = acos(Pos(3,i,j,k)/r);
					double phiAng = atan2(Pos(2,i,j,k),Pos(1,i,j,k));

					// Calculate weighting function

					// Are we in the inner "slope" of Simpson's Rule?
					if (r<simParams.re-simParams.delta.Val(1)/4. && r>simParams.re-5.*simParams.delta.Val(1)/4.)
					{
						w = (r-simParams.re)/simParams.delta.Val(1) + 1.25;
					}
					// Are we in the "flat portion" of Simpson's Rule?
					else if (r>=simParams.re-simParams.delta.Val(1)/4. && r<=simParams.re+simParams.delta.Val(1)/4)
					{
						w = 1.0;
					}
					// Are we in the outer "slope" of Simpson's Rule?
					else if (r>simParams.re+simParams.delta.Val(1)/4. &&  r<simParams.re+5.*simParams.delta.Val(1)/4.)
					{
						w = -(r-simParams.re)/simParams.delta.Val(1) + 1.25;
					}
					// Then we must be outside of the grid shell.
					else
					{
						w = 0.0;
					}
						
					// Using each harmonic array index (which corresponds to a given set of harmonic indices)
					//	generate a rank-2 tensor by calculating the harmonic decoposition of a value of one
					//	for each harmonic
					for (int a=0; a<numHarm(); a++)
					{
						for (int c=0; c<numHarm(); c++)
						{
							// Find the harmonic indices which corresponds to the harmonic array index
							int l1 = getK(0,a);
							int m1 = getK(1,a);
							int n1 = getK(2,a);

							int l2 = getK(0,c);
							int m2 = getK(1,c);
							int n2 = getK(2,c);

							// Calculate the harmonic decomposition using the Misner grid-shell
							//	extraction method
							
							double temp = w*util.misnerF(n1,r,simParams.re,simParams.delta.Val(1)*.75)
								*util.misnerF(n2,r,simParams.re,simParams.delta.Val(1)*.75)
								*simParams.delta.Val(1)*simParams.delta.Val(2)*simParams.delta.Val(3);

							// if m is negative use the imaginary part of the spherical harmonic, otherise
							//	use the real part 
							// CHECK: The imaginary part is still being passed a negative value of m,
							//	is this correct?
							if (m1>=0)
							{
								temp *= real(util.sphericalHarmonicY(l1,m1,theta,phiAng));
							}
							else
							{
								temp *= imag(util.sphericalHarmonicY(l1,m1,theta,phiAng));
							}

							if (m2>=0)
							{
								temp *= real(util.sphericalHarmonicY(l2,m2,theta,phiAng));
							}
							else
							{
								temp *= imag(util.sphericalHarmonicY(l2,m2,theta,phiAng));
							}

							g.Set(g.Val(a,c) + temp,a,c);

						}
					}
				}			// End of central excise section
			}
		}
	}
	return g;
}

void MSField::calcMisnerHarmonics()
{

	// Create a Util object, which has all the various special function (spherical harmonics,
	//	Legendre polynomials, etc.) we will need.
	Util util;
	
	// It's not necessary to scan the grid since the Simpson's Rule weighting falls off
	//	inside and outside of the shell.  I've added a little "Tennessee windage" since the 
	//	Index routine doesn't necessarily give the actual closest grid location to a 
	//	coordinate position.
	int xMin = 0; //Index(1,-(simParams.re+9.*simParams.delta.Val(1)/4.));
	int xMax = (int)simParams.max.Val(1)-1; //Index(1,simParams.re+9.*simParams.delta.Val(1)/4.);
	int yMin = 0; //Index(2,-(simParams.re+9.*simParams.delta.Val(1)/4.));
	int yMax = (int)simParams.max.Val(2)-1; //Index(2,simParams.re+9.*simParams.delta.Val(1)/4.);
	int zMin = 0; //Index(3,-(simParams.re+9.*simParams.delta.Val(1)/4.));
	int zMax = (int)simParams.max.Val(3)-1; //Index(3,simParams.re+9.*simParams.delta.Val(1)/4.);
		
	// Initialize the harmonic coefficient array	
	for (int a = 0; a<numHarm(); a++)
	{
		b[a] = 0.;
	}
	
	// Scan within the shell
	for (int i = xMin; i<=xMax; i++)
	{
		for (int j = yMin; j<=yMax; j++)
		{
			for (int k = zMin; k<=zMax; k++)
			{
				double w = 0.0;

				// Calculate the coordinates of the grid location in spherical coordinates
				double r = rDist(i,j,k);
				if (r>=simParams.re-5.*simParams.delta.Val(1)/4. && r<=simParams.re+5.*simParams.delta.Val(1)/4.)
				{
					double theta = acos(Pos(3,i,j,k)/r);
					double phiAng = atan2(Pos(2,i,j,k),Pos(1,i,j,k));

					// Calculate weighting function
					// Are we in the inner "slope" of Simpson's Rule?
					if (r<simParams.re-simParams.delta.Val(1)/4. && r>simParams.re-5.*simParams.delta.Val(1)/4.)
					{
						w = (r-simParams.re)/simParams.delta.Val(1) + 1.25;
					}
					// Are we in the "flat portion" of Simpson's Rule?
					else if (r>=simParams.re-simParams.delta.Val(1)/4. && r<=simParams.re+simParams.delta.Val(1)/4)
					{
						w = 1.0;
					}
					// Are we in the outer "slope" of Simpson's Rule?
					else if (r>simParams.re+simParams.delta.Val(1)/4. &&  r<simParams.re+5.*simParams.delta.Val(1)/4.)
					{
						w = -(r-simParams.re)/simParams.delta.Val(1) + 1.25;
					}
					// Then we must be outside of the grid shell.
					else
					{
						w = 0.0;
					}

					// For each harmonic array index, calculate the harmonic cooefficient.
					//	We must also contract with the normalization tensor, so we sum as well. 
					for (int a=0; a<numHarm(); a++)
					{
						for (int c=0; c<numHarm(); c++)
						{		
							// Calculate the harmonic indices corresponding to the inner loop's
							//	harmonic array index
							int l = getK(0,c);
							int m = getK(1,c);
							int n = getK(2,c);
							
							// Calculate the harmonic decomposition using the Misner grid-shell
							//	extraction method, normalized by the normalization tensor.
							// If m is negative use the imaginary part of the spherical harmonic, otherise
							//	use the real part 
							// CHECK: The imaginary part is still being passed a negative value of m,
							//	is this correct?
							if (m>=0)
							{
								b[a] = b[a] + w*g.Val(a,c)*phi[1][i][j][k]
									*real(util.sphericalHarmonicY(l,m,theta,phiAng))
									*util.misnerF(n,r,simParams.re,simParams.delta.Val(1)*0.75)
									*simParams.delta.Val(1)*simParams.delta.Val(2)*simParams.delta.Val(3);
							}
							else
							{
								b[a] = b[a] + w*g.Val(a,c)*phi[1][i][j][k]
									*imag(util.sphericalHarmonicY(l,m,theta,phiAng))
									*util.misnerF(n,r,simParams.re,simParams.delta.Val(1)*0.75)
									*simParams.delta.Val(1)*simParams.delta.Val(2)*simParams.delta.Val(3);
							}
						}
					}
				}			// End of central excise condition
			}
		}
	}
	
	return;
}


// At a given position, reconstruct phi from the harmonic decomposition
// NOTE: This routine, applied over the grid, slows the code down considerably.
//	Normally, one is calculating phi at a grid location so the alternative form
//	of this routine (see below) should be used.  However, in order to reconstruct
//	the field off of the grid (e.g., at LIGO), this is the routine which must be 
//	used.
double MSField::restoreFromMisnerHarmonics(tensor pos)
{
	double output = 0.0;
	int a;
		
	// Create a Util object, which has all the various special function (spherical harmonics,
	//	Legendre polynomials, Bessel functions, etc.) we will need.
	Util util;

	// Calculate the coordinates of the requested position in spherical coordinates
	double r = pos.Vmag();
	double theta = 0.;
	double phiAng = 0.;
	if (r!=0.)
	{
		theta = acos(pos.Val(3)/r);
		phiAng = atan2(pos.Val(2),pos.Val(1));
	}
	else
	{
		cout << "WARNING: Central coordinate singularity was included in harmonic restore. \n";
	}

		
	// Calculate the sums we need
	array gridSum(1,numHarm());
		
	// Multiply the harmonic coefficient by the sum of the Misner harmonic terms we are using
	//	Note that the data structures are set up so that all the misner harmonics are 
	//	sequential in order to facilitate this sum
	for (a=0; a<numHarm(); a=a+simParams.maxN+1)
	{
		gridSum.Set(b[a]*util.misnerF(0,simParams.re,simParams.re,simParams.delta.Val(1)*0.75)
			+ b[a+1]*util.misnerF(1,simParams.re,simParams.re,simParams.delta.Val(1)*0.75)
			+ b[a+2]*util.misnerF(2,simParams.re,simParams.re,simParams.delta.Val(1)*0.75)
			,a);
	}
				
	// Now reconstruct from the spherical harmonics
	//	Note that we have already summed over the Misner harmonics and only that sum is needed here.
	for (a=0; a<numHarm(); a=a+simParams.maxN+1)
	{
		// Get the harmonic indices which correspond to the harmonic array index
		int l = getK(0,a);
		int m = getK(1,a);
				
		// Reconstruct from the analytic solution to the field equation (valid where the stress-energy
		//	tensor is zero).  Remember that the m=0 terms have to be handled seperately from the others.
		if (m==0)
		{
			output = output + pow(simParams.re/r,l+1)*real(util.sphericalHarmonicY(l,0,theta,phiAng))
				*gridSum.Val(a);
		}
		// If m is negative use the imaginary part of the spherical harmonic, otherise
		//	use the real part 
		// CHECK: The imaginary part is still being passed a negative value of m,
		//	is this correct?
		else if (m>0)
		{
			// I'm assuming that for positive m we take the real part, as this was what was 
			//	done before adding the sphBesJ - be sure to verify this
			complex tempOutput;
			tempOutput = (complex(util.sphNeuN(l,(double)m*simParams.omega*r),
				-util.sphBesJ(l,(double)m*simParams.omega*r))) /
				(complex(util.sphNeuN(l,(double)m*simParams.omega*simParams.re),
				-util.sphBesJ(l,(double)m*simParams.omega*simParams.re)));

			output = output + real(tempOutput*util.sphericalHarmonicY(l,m,theta,phiAng))
				*gridSum.Val(a);
		}
		else
		{
			// I'm assuming that for negative m we take the imag part, as this was what was 
			//	done before adding the sphBesJ - be sure to verify this
			complex tempOutput;
			tempOutput = (complex(util.sphNeuN(l,(double)m*simParams.omega*r),
				-util.sphBesJ(l,(double)m*simParams.omega*r))) /
				(complex(util.sphNeuN(l,(double)m*simParams.omega*simParams.re),
				-util.sphBesJ(l,(double)m*simParams.omega*simParams.re)));

			output = output + imag(tempOutput*util.sphericalHarmonicY(l,m,theta,phiAng))
				*gridSum.Val(a);
		}
	}
					
	return output;
}

// At a given position, reconstruct phi from the harmonic decomposition
double MSField::restoreFromMisnerHarmonics(int i, int j, int k)
{
	double output = 0.0;
	int a;

	// Create a Util object, which has all the various special function (spherical harmonics,
	//	Legendre polynomials, Bessel functions, etc.) we will need.
	Util util;
		
	// Calculate the coordinates of the requested position in spherical coordinates
	double r = rDist(i,j,k);
	double theta = 0.;
	double phiAng = 0.;
	if (r!=0.)
	{
		theta = acos(Pos(3,i,j,k)/r);
		phiAng = atan2(Pos(2,i,j,k),Pos(1,i,j,k));
	}
	else
	{
		cout << "WARNING: Central coordinate singularity was included in harmonic restore. \n";
	}
		
	// Calculate the sums we need
	array gridSum(1,numHarm());
		
	// Multiply the harmonic coefficient by the sum of the Misner harmonic terms we are using
	//	Note that the data structures are set up so that all the misner harmonics are 
	//	sequential in order to facilitate this sum
	for (a=0; a<numHarm(); a=a+simParams.maxN+1)
	{
		gridSum.Set(b[a]*util.misnerF(0,simParams.re,simParams.re,simParams.delta.Val(1)*0.75)
			+ b[a+1]*util.misnerF(1,simParams.re,simParams.re,simParams.delta.Val(1)*0.75)
			+ b[a+2]*util.misnerF(2,simParams.re,simParams.re,simParams.delta.Val(1)*0.75)
			,a);
	}

	// Now reconstruct from the spherical harmonics
	//	Note that we have already summed over the Misner harmonics and only that sum is needed here.
	for (a=0; a<numHarm(); a=a+simParams.maxN+1)
	{
		// Get the harmonic indices which correspond to the harmonic array index
		int l = getK(0,a);
		int m = getK(1,a);
				
		// Reconstruct from the analytic solution to the field equation (valid where the stress-energy
		//	tensor is zero).  Remember that the m=0 terms have to be handled seperately from the others.
		if (m==0)
		{
			output = output + pow(simParams.re/r,l+1)*real(util.sphericalHarmonicY(l,0,theta,phiAng))
				*gridSum.Val(a);
		}
		// If m is negative use the imaginary part of the spherical harmonic, otherise
		//	use the real part 
		// CHECK: The imaginary part is still being passed a negative value of m,
		//	is this correct?
		else if (m>0)
		{
			// I'm assuming that for positive m we take the real part, as this was what was 
			//	done before adding the sphBesJ - be sure to verify this
			complex tempOutput;
			tempOutput = (complex(util.sphNeuN(l,(double)m*simParams.omega*r),
				-util.sphBesJ(l,(double)m*simParams.omega*r))) /
				(complex(util.sphNeuN(l,(double)m*simParams.omega*simParams.re),
				-util.sphBesJ(l,(double)m*simParams.omega*simParams.re)));

			output = output + real(tempOutput*util.sphericalHarmonicY(l,m,theta,phiAng))
				*gridSum.Val(a);
		}
		else
		{
			// I'm assuming that for negative m we take the imag part, as this was what was 
			//	done before adding the sphBesJ - be sure to verify this
			complex tempOutput;
			tempOutput = (complex(util.sphNeuN(l,(double)m*simParams.omega*r),
				-util.sphBesJ(l,(double)m*simParams.omega*r))) /
				(complex(util.sphNeuN(l,(double)m*simParams.omega*simParams.re),
				-util.sphBesJ(l,(double)m*simParams.omega*simParams.re)));

			output = output + imag(tempOutput*util.sphericalHarmonicY(l,m,theta,phiAng))
				*gridSum.Val(a);
		}
	}
					
	return output;
}



// This routine simply combines all the harmonic decomposition intialization
//	tasks in one place

void MSField::initMisnerHarmonics()
{
	cout << "Initializing spherical harmonics...\n";
		
	// Create a Util object, which has the invert function we will need.
	Util util;

	int i,j;
		
	// Resize the harmonic coefficents array to the correct number of harmonics
	cout << "Resizing harmonic array...\n";
	b = new double [numHarm()];

	// In the Java version we load the normalization tensor from an ASCII 
	//	file called "g.dat".  If the first value is "-1", the routine calculates
	//	a new normalization tensor, otherwise it just loads the stored one.
	//	I'll that capability to the C++ version eventually.
	cout << "Calculating new normalization tensor...\n";
	g = calcNormTensor();			

	// Invert the normalization tensor so that we don't have to define tensor division(!)
	cout << "Inverting matrix...\n";
	util.invert(g);
		
	// write the inverted normalization tensor to an output file called "gInv.dat" 
	//	This file can be renamed to "g.dat" and used in the manner described above.
	FILE *gFile;
	gFile = fopen("gInv.dat","w");

	for (i=0;i<numHarm();i++)
	{
		for (j=0;j<numHarm();j++)
		{
			fprintf(gFile,"%f",g.Val(i,j));
			fprintf(gFile,"\n");
		}
	}

	// Calculate the spherical harmonics decomposition
	cout << "Calculating initial harmonic decomposition...\n";
	calcMisnerHarmonics();
	cout << "Harmonic initialization complete.\n";
	
}

double MSField::calcDipolePoyntingFlux()
{
	double a1 = 0.;
	double a3 = 0.;
	Util util;

	for (int n=0; n<simParams.maxN; n++)
	{
		a1 += b[getK(1,-1,n)] * util.misnerF(n,simParams.re,simParams.re,simParams.delta.Val(1)*0.75);
		a3 += b[getK(1,+1,n)] * util.misnerF(n,simParams.re,simParams.re,simParams.delta.Val(1)*0.75);
	}


	return pow(C,5)/(4.*PI*G) * pow(simParams.omega*simParams.re,2) * (a1*a3);
}

double MSField::calcQuadPoyntingFlux()
{
	double a1 = 0.;
	double a2 = 0.;
	double a4 = 0.;
	double a5 = 0.;
	Util util;

	for (int n=0; n<simParams.maxN; n++)
	{
		a1 += b[getK(2,-2,n)] * util.misnerF(n,simParams.re,simParams.re,simParams.delta.Val(1)*0.75);
		a2 += b[getK(2,-1,n)] * util.misnerF(n,simParams.re,simParams.re,simParams.delta.Val(1)*0.75);
		a4 += b[getK(2,+1,n)] * util.misnerF(n,simParams.re,simParams.re,simParams.delta.Val(1)*0.75);
		a5 += b[getK(2,+2,n)] * util.misnerF(n,simParams.re,simParams.re,simParams.delta.Val(1)*0.75);
	}


	return pow(C,5)/(4.*PI*G) * pow(simParams.omega*simParams.re,2) * (a1*a5 - 4.*a2*a4);
}

#ifdef PART
void MSField::RegisterParticles(Particle *particle)
{
 	my_particle_partner = particle;
}
#endif

// (KW) These two functions are never called, so I've commented them out completely.
//			They were used to generate Charlie's test data for his Santa Barbara talk.
/*double MSField::matterVar(int q, int x, int y, int z, int c, int d)
{

	// q=timeslice (old, current, or new)
	// x,y,z=grid location
	// c,d=requested component

	double dLdg = 0.;
/*
	// Add up the contribution from each particle
	for (int a=0; a<2; a++)
	{
		// Calculate the squared distance from the particle's center
		//	to the grid location specified as a fraction of the
		//	smoothing length squared
		double r2 = (pow(particle.Position(a,1)-Pos(1,x,y,z),2.)
				+ pow(particle.Position(a,2)-Pos(2,x,y,z),2.)
				+ pow(particle.Position(a,3)-Pos(3,x,y,z),2.))
			/pow(particle.SmoothLength(a,1),2.0);

		// w is the weighting kernel
		double w = 0.0;

		// norm is the normalization of the weighting kernel
		double norm = 1.0;
		if (r2 < 1.0)
		{
			// Calculate normalization constant	by "smoothing" a value of 1
			norm = 0.0;
			double rnorm2;

			// It's not necessary to scan the entire grid, since
			//	anything outside the smoothing length is zero.
			// NOTE: It would probably be a good idea to add some
			//	"Tennessee windage" to this since the Index routine
			//	may not reliably giving the closest grid location.
			int minX = Index(1,particle.Position(a,1)
				-particle.SmoothLength(a,1)*1.5);
			int maxX = Index(1,particle.Position(a,1)
				+particle.SmoothLength(a,1)*1.5);
			int minY = Index(2,particle.Position(a,2)
				-particle.SmoothLength(a,2)*1.5);
			int maxY = Index(2,particle.Position(a,2)
				+particle.SmoothLength(a,2)*1.5);
			int minZ = Index(3,particle.Position(a,3)
				-particle.SmoothLength(a,3)*1.5);
			int maxZ = Index(3,particle.Position(a,3)
				+particle.SmoothLength(a,3)*1.5);

			// Scan within the mins and maxs defined above
			for (int i=minX; i<=maxX; i++)
			{
				for (int j=minY; j<=maxY; j++)
				{
					for (int k=minZ; k<=maxZ; k++)
					{
						// Calculate the squared distance between the particle
						//	center and the grid location being scanned as a
						//	fraction of the smoothing length squared

						rnorm2 = (pow(particle.Position(a,1)-Pos(1,i,j,k),2.)
								+ pow(particle.Position(a,2)-Pos(2,i,j,k),2.)
								+ pow(particle.Position(a,3)-Pos(3,i,j,k),2.))
							/pow(particle.SmoothLength(a,1),2.0);

						if (rnorm2 < 1.0)
						{
							// The smoothing kernel Charlie and I worked out
							//	Note that we are still "smoothing" a value of 1 here
							norm  = norm + pow(1.0-rnorm2,simParams.exponent);
						}
					}
				}
			}
			// Calculate the actual weighting of the requested grid location
			w = pow(1.0-r2,simParams.exponent);
		}

		// Calculate the contribution of this particle to the total density at the
		//	requested location using Conrad's density definition
		dLdg += w/norm*particle.Mass(a)/2.
			* metricG(q,x,y,z,0,c)*metricG(q,x,y,z,0,d)/sqrt(-metricG(q,x,y,z,0,0));

	}

	return dLdg;
}

// Calculate the total smoothed density at any grid location
//	usuing Conrad's SPH
double MSField::calcDensity(int q, int x, int y, int z)
{

	// q=timeslice (old, current, or new)
	// x,y,z=grid location

	double totalRho = 0.;
/*
	// Add up the contribution from each particle
	for (int a=0; a<2; a++)
	{
		// Calculate the squared distance from the particle's center
		//	to the grid location specified as a fraction of the
		//	smoothing length squared
		double r2 = (pow(particle.Position(a,1)-Pos(1,x,y,z),2.)
				+ pow(particle.Position(a,2)-Pos(2,x,y,z),2.)
				+ pow(particle.Position(a,3)-Pos(3,x,y,z),2.))
			/pow(particle.SmoothLength(a,1),2.0);

		// w is the weighting kernel
		double w = 0.0;

		// norm is the normalization of the weighting kernel
		double norm = 1.0;
		if (r2 < 1.0)
		{
			// Calculate normalization constant	by "smoothing" a value of 1
			norm = 0.0;
			double rnorm2;

			// It's not necessary to scan the entire grid, since
			//	anything outside the smoothing length is zero.
			// NOTE: It would probably be a good idea to add some
			//	"Tennessee windage" to this since the Index routine
			//	may not reliably giving the closest grid location.
			int minX = Index(1,particle.Position(a,1)
				-particle.SmoothLength(a,1)*1.5);
			int maxX = Index(1,particle.Position(a,1)
				+particle.SmoothLength(a,1)*1.5);
			int minY = Index(2,particle.Position(a,2)
				-particle.SmoothLength(a,2)*1.5);
			int maxY = Index(2,particle.Position(a,2)
				+particle.SmoothLength(a,2)*1.5);
			int minZ = Index(3,particle.Position(a,3)
				-particle.SmoothLength(a,3)*1.5);
			int maxZ = Index(3,particle.Position(a,3)
				+particle.SmoothLength(a,3)*1.5);

			// Scan within the mins and maxs defined above
			for (int i=minX; i<=maxX; i++)
			{
				for (int j=minY; j<=maxY; j++)
				{
					for (int k=minZ; k<=maxZ; k++)
					{
						// Calculate the squared distance between the particle
						//	center and the grid location being scanned as a
						//	fraction of the smoothing length squared

						rnorm2 = (pow(particle.Position(a,1)-Pos(1,i,j,k),2.)
								+ pow(particle.Position(a,2)-Pos(2,i,j,k),2.)
								+ pow(particle.Position(a,3)-Pos(3,i,j,k),2.))
							/pow(particle.SmoothLength(a,1),2.0);

						if (rnorm2 < 1.0)
						{
							// The smoothing kernel Charlie and I worked out
							//	Note that we are still "smoothing" a value of 1 here
							norm  = norm + pow(1.0-rnorm2,simParams.exponent);
						}
					}
				}
			}
			// Calculate the actual weighting of the requested grid location
			w = pow(1.0-r2,simParams.exponent);
		}

		// Calculate the contribution of this particle to the total density at the
		//	requested location using a smoothed density (Note: This formula is not
		//  even close to being correct, it is strictly for testing the scalar code.
		totalRho += w/norm*particle.Mass(a)
//			* metricG(q,x,y,z,0,0)*metricG(q,x,y,z,0,0)/sqrt(-metricG(q,x,y,z,0,0));
	        * sqrt(exp(6.0*phi.Val(q,x,y,z)));

	}

	return totalRho;
}
*/
