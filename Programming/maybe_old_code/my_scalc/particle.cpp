// Particle.cpp
// Encapsulates a set of point particles - this object will contain all of Conrad's code
// Keith Watt 1999

#include "Particle.hpp"

//******************************************************************************
//Name:  Particle                                                              *
//                                                                             *
//Purpose:  the particle container class constructor (default)                 *
//                                                                             *
//Takes:                                                                       *
//******************************************************************************
Particle::Particle()
{
	numParticles = 2;
   for(int rk_stage = 0; rk_stage < 3; rk_stage++)
    {
      //position and velocities
      pos[rk_stage]              = new tensor [numParticles];
   	u_dn[rk_stage]             = new tensor [numParticles];
      u_up[rk_stage]             = new tensor [numParticles];
      u0_dn[rk_stage]            = new double [numParticles];
      u0_up[rk_stage]            = new double [numParticles];
      four_U_up[rk_stage]        = new tensor [numParticles];
      four_U_dn[rk_stage]        = new tensor [numParticles];
      //inverse 4-metric and its derivative
      inv_g[rk_stage]            = new tensor [numParticles];
      d_inv_g[rk_stage]          = new tensor [numParticles];
      //ADM decomposition
      lapse[rk_stage]            = new double [numParticles];
      shift[rk_stage]            = new tensor [numParticles];
      three_metric[rk_stage]     = new tensor [numParticles];
      inv_three_metric[rk_stage] = new tensor [numParticles];
      norm[rk_stage]             = new double [numParticles];
      d_norm[rk_stage]           = new tensor [numParticles];
    }
	masses            = new double [numParticles];
   smoothing_lengths = new double [numParticles];
   T_temp4x4x3.Resize(3,4,4,3);
   T_temp4x4.Resize(2,4,4);
   T_temp3x3.Resize(2,3,3);

   for(int i = 0; i < numParticles; i++)
    {
     for(int rk_stage = 0; rk_stage < 3; rk_stage++)   //loop over the RK stages
      {
       sprintf(name,"pos%d_%d",rk_stage,i);
       pos[rk_stage][i].SetName(name);
       sprintf(name,"udn%d_%d",rk_stage,i);
       u_dn[rk_stage][i].SetName(name);
       sprintf(name,"uup%d_%d",rk_stage,i);
       u_up[rk_stage][i].SetName(name);
       sprintf(name,"4up%d_%d",rk_stage,i);
       four_U_up[rk_stage][i].SetName(name);
       sprintf(name,"4dn%d_%d",rk_stage,i);
       four_U_dn[rk_stage][i].SetName(name);
       sprintf(name,"ig%d_%d",rk_stage,i);
       inv_g[rk_stage][i].SetName(name);
       sprintf(name,"dig%d_%d",rk_stage,i);
       d_inv_g[rk_stage][i].SetName(name);
       sprintf(name,"shf%d_%d",rk_stage,i);
       shift[rk_stage][i].SetName(name);
       sprintf(name,"3g%d_%d",rk_stage,i);
       three_metric[rk_stage][i].SetName(name);
       sprintf(name,"i3g%d_%d",rk_stage,i);
       inv_three_metric[rk_stage][i].SetName(name);
       sprintf(name,"dW%d_%d",rk_stage,i);
       d_norm[rk_stage][i].SetName(name);
       pos[rk_stage][i].Resize(1,3);
       u_dn[rk_stage][i].Resize(1,3);
       u_up[rk_stage][i].Resize(1,3);
       four_U_up[rk_stage][i].Resize(1,4);
       four_U_dn[rk_stage][i].Resize(1,4);
	    inv_g[rk_stage][i].Resize(2,4,4);
	    d_inv_g[rk_stage][i].Resize(3,4,4,3);
       shift[rk_stage][i].Resize(1,3);
       three_metric[rk_stage][i].Resize(2,3,3);
       inv_three_metric[rk_stage][i].Resize(2,3,3);
       d_norm[rk_stage][i].Resize(1,3);
      }
    }
  last_update_time = -1.0;
  last_stage       = -1;
  Delta_mag = 0.0;
  y = 0.0;
  for(int i = 0; i < 3; i++)
  {
	  Deltas[i]           = 0.0;
	  current_pos[i]      = 0.0;
	  reference_pos[i]    = 0.0;
     pos_diff[i]         = 0.0;
     spatial_extent[i]   = 0.0;
     spatial_extent[i+3] = 0.0;
  }

  particle_smoothing_method = 0;
  grid_smoothing_method = 0;
  out_of_range = FALSE;
  is_periodic = FALSE;
  grid_radius = 0;
  grid_upper_bound = 0;
  grid_lower_bound = 0;
  smoothing_length = 1.0;
  normalization = 0.0;
  is_series_set = FALSE;
  is_deriv_field = FALSE;
  my_field_partner = NULL;

  return;
}

//******************************************************************************
//Name:  Particle                                                              *
//                                                                             *
//Purpose:  the particle container class constructor                           *
//                                                                             *
//Takes:    the number of particles                                            *
//******************************************************************************
Particle::Particle(int number_of_particles)
{
   for(int rk_stage = 0; rk_stage < 3; rk_stage++)
    {
      //position and velocities
      pos[rk_stage]              = new tensor [numParticles];
   	u_dn[rk_stage]             = new tensor [numParticles];
      u_up[rk_stage]             = new tensor [numParticles];
      u0_dn[rk_stage]            = new double [numParticles];
      u0_up[rk_stage]            = new double [numParticles];
      four_U_up[rk_stage]        = new tensor [numParticles];
      four_U_dn[rk_stage]        = new tensor [numParticles];
      //inverse 4-metric and its derivative
      inv_g[rk_stage]            = new tensor [numParticles];
      d_inv_g[rk_stage]          = new tensor [numParticles];
      //ADM decomposition
      lapse[rk_stage]            = new double [numParticles];
      shift[rk_stage]            = new tensor [numParticles];
      three_metric[rk_stage]     = new tensor [numParticles];
      inv_three_metric[rk_stage] = new tensor [numParticles];
      norm[rk_stage]             = new double [numParticles];
      d_norm[rk_stage]           = new tensor [numParticles];
    }
	masses            = new double [numParticles];
   smoothing_lengths = new double [numParticles];
   T_temp4x4x3.Resize(3,4,4,3);
   T_temp4x4.Resize(2,4,4);
   T_temp3x3.Resize(2,3,3);

   for(int i = 0; i < numParticles; i++)
    {
     for(int rk_stage = 0; rk_stage < 3; rk_stage++)   //loop over the RK stages
      {
       sprintf(name,"pos%d_%d",rk_stage,i);
       pos[rk_stage][i].SetName(name);
       sprintf(name,"udn%d_%d",rk_stage,i);
       u_dn[rk_stage][i].SetName(name);
       sprintf(name,"uup%d_%d",rk_stage,i);
       u_up[rk_stage][i].SetName(name);
       sprintf(name,"4up%d_%d",rk_stage,i);
       four_U_up[rk_stage][i].SetName(name);
       sprintf(name,"4dn%d_%d",rk_stage,i);
       four_U_dn[rk_stage][i].SetName(name);
       sprintf(name,"ig%d_%d",rk_stage,i);
       inv_g[rk_stage][i].SetName(name);
       sprintf(name,"dig%d_%d",rk_stage,i);
       d_inv_g[rk_stage][i].SetName(name);
       sprintf(name,"shf%d_%d",rk_stage,i);
       shift[rk_stage][i].SetName(name);
       sprintf(name,"3g%d_%d",rk_stage,i);
       three_metric[rk_stage][i].SetName(name);
       sprintf(name,"i3g%d_%d",rk_stage,i);
       inv_three_metric[rk_stage][i].SetName(name);
       sprintf(name,"dW%d_%d",rk_stage,i);
       d_norm[rk_stage][i].SetName(name);
       pos[rk_stage][i].Resize(1,3);
       u_dn[rk_stage][i].Resize(1,3);
       u_up[rk_stage][i].Resize(1,3);
       four_U_up[rk_stage][i].Resize(1,4);
       four_U_dn[rk_stage][i].Resize(1,4);
	    inv_g[rk_stage][i].Resize(2,4,4);
	    d_inv_g[rk_stage][i].Resize(3,4,4,3);
       shift[rk_stage][i].Resize(1,3);
       three_metric[rk_stage][i].Resize(2,3,3);
       inv_three_metric[rk_stage][i].Resize(2,3,3);
       d_norm[rk_stage][i].Resize(1,3);
      }
    }
  last_update_time = -1.0;
  last_stage       = -1;
  particle_smoothing_method = 0;
  grid_smoothing_method = 0;
  out_of_range = FALSE;
  is_periodic = FALSE;
  grid_radius = 0;
  grid_upper_bound = 0;
  grid_lower_bound = 0;
  smoothing_length = 1.0;
  normalization = 0.0;
  is_series_set = FALSE;
  is_deriv_field = FALSE;
  my_field_partner = NULL;

  return;
}

//******************************************************************************
//Name:  RegisterField                                                         *
//                                                                             *
//Purpose:  to make the location of the field known to the particle class      *
//                                                                             *
//Takes: a pointer to the field                                                *
//******************************************************************************
void Particle::RegisterField(MSField *field)
{
 	my_field_partner = field;
}


//******************************************************************************
//Name:  matterVar                                                             *
//                                                                             *
//Purpose:  returns the value of \frac{ \partial L}{\partial g^{ab} }          *
//                                                                             *
//Takes: an int telling what stage of the RK and 3 ints telling the grid       *
//       location                                                              *
//******************************************************************************
tensor Particle::matterVar(int stage, int i, int j, int k)
{

	//update the particle data if first time called
	if( stage != last_stage )
	{
		update_inv_g(stage);
      ADM_decomposition(stage);
      update_u_dns(stage);
      update_u_ups(stage);
      RK_step(stage);
      last_stage = stage;
   }

	//now return T_{\mu\nu}
   indices[0] = i;
   indices[1] = j;
   indices[2] = k;
   grid_pos <= GetPosition(indices);
   T_temp4x4 <= 0.0*T_temp4x4;
	for ( c1 = 0; c1 < numParticles; c1++)
	{
		pos_diff_tensor <= grid_pos - pos[stage][c1];
      y = pos_diff_tensor.Vmag() / smoothing_length;
      W = Smoother(CWM_SMOOTHER);
      //only calculate the matter variation term if within smoothing radius
      if( fabs(W) > 0 + SMOOTHING_TOL)
       {
 		   factor  = -0.5 * masses[c1] / u0_up[stage][c1];
        	T_temp4x4 <= T_temp4x4 + factor * four_U_dn[stage][c1]
                                         * four_U_dn[stage][c1]
                                         * W;
       }
	}

    //Currently stubbed out to return the identity matrix (more or less)
    return T_temp4x4;

}


//******************************************************************************
//Name:  update_inv_g                                                          *
//                                                                             *
//Purpose:  scans over the relevant regions of the grid and gets a smoothed    *
//          version of the inverse metric at the particle                      *
//                                                                             *
//Takes:                                                                       *
//******************************************************************************
void Particle::update_inv_g(int stage)
{
 int i, j, k, m, n;
 double scaled_pos[3];
 int grid_pos[3];
 double position[3];
 double W;

//loop over the particles and determine the smoothed inverse metric
//for each one
 for( c1 = 0; c1 < numParticles; c1++)
 {
  //Set the reference position using the particle's position
  reference_pos[0] = pos[stage][c1].Val(0);
  reference_pos[1] = pos[stage][c1].Val(1);
  reference_pos[2] = pos[stage][c1].Val(2);

  //Get the real indices for the scaled position
  GetIndices(reference_pos, scaled_pos);

  //Reset the normalization
  normalization = 0.0;

  //initialize member data tensors
  sum_Wm      <= 0.0*sum_Wm;
  T_temp4x4   <= 0.0*T_temp4x4;
  T_temp4x4x3 <= 0.0*T_temp4x4x3;

  //step through the floor steping
  for( i = grid_lower_bound; i < grid_upper_bound; i++)
   {
	 for( j = grid_lower_bound; j < grid_upper_bound; j++)
	  {
	   for( k = grid_lower_bound; k < grid_upper_bound; k++)
		 {
         //Get the indices of the lowest corner of the cell
			grid_pos[0] = (int)floor( scaled_pos[0] + i );
			grid_pos[1] = (int)floor( scaled_pos[1] + j );
         grid_pos[2] = (int)floor( scaled_pos[2] + k );

         //find the metric at the particular grid point
         inv_g_temp <= my_field_partner->GetInvG(stage,
                                                 grid_pos[0],
                                                 grid_pos[1],
                                                 grid_pos[2]);

         //find the position this grid point
         GetPosition(grid_pos,current_pos);

         //form the difference array and its magnitude (for use by smoother)
         //and the difference position
         y = 0;
         for( m = 0; m < 3; m++ )
			 {
           pos_diff[m] = reference_pos[m] - current_pos[m];
           y += pos_diff[m]*pos_diff[m];
			 }
         y = sqrt(y)/smoothing_length;
         if( y > 1.0 + SMOOTHING_TOL ) continue;

         //get the (psuedo) scalar smoothing value
         W = Smoother(CWM_SMOOTHER);

         //get the derivative of the smoothing
         Smoother_Deriv(CWM_SMOOTHER);

         //keep running sum
         sum_Wm      <= sum_Wm + Wm;
 	      T_temp4x4   <= T_temp4x4 + W*inv_g_temp;
         T_temp4x4x3 <= T_temp4x4x3 + inv_g_temp*Wm;

	    }
	  }
   }

  //Normalize the interpolated tensor and pack arrays
  inv_g[stage][c1]   <= T_temp4x4;
  d_inv_g[stage][c1] <= T_temp4x4x3;
  norm[stage][c1]     = normalization;
  d_norm[stage][c1]  <= sum_Wm;

 }



}

//******************************************************************************
//Name:  ADM_decomposition                                                     *
//                                                                             *
//Purpose:  returns the value of T_{\mu\nu} at a specified grid point          *
//                                                                             *
//Takes: an int telling what stage of the RK and 3 ints telling the grid       *
//       location                                                              *
//******************************************************************************
void Particle::ADM_decomposition(int stage)
{
  for(c1 = 0; c1 < numParticles; c1++)
   {
     //take the 0,0 component of g_inv to get the lapse
     lapse[stage][c1] = sqrt(-1.0/inv_g[stage][c1].p->m[0]);
     inv_lapse_val = 1.0/lapse[stage][c1];
     inv_lapse_sqr = inv_lapse_val*inv_lapse_val;

     //take the spatial portion of the top row of inv_g for the shift
     shift[stage][c1].p->m[0] = inv_g[stage][c1].p->m[1]*inv_lapse_sqr;
     shift[stage][c1].p->m[1] = inv_g[stage][c1].p->m[2]*inv_lapse_sqr;
     shift[stage][c1].p->m[2] = inv_g[stage][c1].p->m[3]*inv_lapse_sqr;

     //take the spatial portion of g_inv for the inverse 3-metric
     T_temp3x3 <= inv_lapse_sqr*shift[stage][c1]*shift[stage][c1];

     inv_three_metric[stage][c1].p->m[0] =   inv_g[stage][c1].p->m[5]
                                           + T_temp3x3.p->m[0];
     inv_three_metric[stage][c1].p->m[1] =   inv_g[stage][c1].p->m[6]
                                           + T_temp3x3.p->m[1];
     inv_three_metric[stage][c1].p->m[2] =   inv_g[stage][c1].p->m[7]
                                           + T_temp3x3.p->m[2];
     inv_three_metric[stage][c1].p->m[3] =   inv_g[stage][c1].p->m[9]
                                           + T_temp3x3.p->m[3];
     inv_three_metric[stage][c1].p->m[4] =   inv_g[stage][c1].p->m[10]
                                           + T_temp3x3.p->m[4];
     inv_three_metric[stage][c1].p->m[5] =   inv_g[stage][c1].p->m[11]
                                           + T_temp3x3.p->m[5];
     inv_three_metric[stage][c1].p->m[6] =   inv_g[stage][c1].p->m[13]
                                           + T_temp3x3.p->m[6];
     inv_three_metric[stage][c1].p->m[7] =   inv_g[stage][c1].p->m[14]
                                           + T_temp3x3.p->m[7];
     inv_three_metric[stage][c1].p->m[8] =   inv_g[stage][c1].p->m[15]
                                           + T_temp3x3.p->m[8];

     //get the determinant of the inverse three metric prior to forming
     //the inverse
     det = inv_three_metric[stage][c1].p->m[0] *
           ( inv_three_metric[stage][c1].p->m[4] *
             inv_three_metric[stage][c1].p->m[8] -
             inv_three_metric[stage][c1].p->m[5] *
             inv_three_metric[stage][c1].p->m[7] )
           +
           inv_three_metric[stage][c1].p->m[1] *
           ( inv_three_metric[stage][c1].p->m[5] *
             inv_three_metric[stage][c1].p->m[6] -
             inv_three_metric[stage][c1].p->m[3] *
             inv_three_metric[stage][c1].p->m[8] )
           +
           inv_three_metric[stage][c1].p->m[2] *
           ( inv_three_metric[stage][c1].p->m[3] *
             inv_three_metric[stage][c1].p->m[7] -
             inv_three_metric[stage][c1].p->m[4] *
             inv_three_metric[stage][c1].p->m[6] );

    //now form the components of the three metric -- this is long and
    //CS's notes should be consulted

    //0,0 Component
    three_metric[stage][c1].p->m[0] = ( inv_three_metric[stage][c1].p->m[4] *
                                        inv_three_metric[stage][c1].p->m[8] -
                                        inv_three_metric[stage][c1].p->m[5] *
                                        inv_three_metric[stage][c1].p->m[7] );
    //0,1 Component
    three_metric[stage][c1].p->m[1] = ( inv_three_metric[stage][c1].p->m[2] *
                                        inv_three_metric[stage][c1].p->m[7] -
                                        inv_three_metric[stage][c1].p->m[1] *
                                        inv_three_metric[stage][c1].p->m[8] );
    //0,2 Component
    three_metric[stage][c1].p->m[2] = ( inv_three_metric[stage][c1].p->m[1] *
                                        inv_three_metric[stage][c1].p->m[5] -
                                        inv_three_metric[stage][c1].p->m[2] *
                                        inv_three_metric[stage][c1].p->m[4] );
    //1,0 Component
    three_metric[stage][c1].p->m[3] = ( inv_three_metric[stage][c1].p->m[5] *
                                        inv_three_metric[stage][c1].p->m[6] -
                                        inv_three_metric[stage][c1].p->m[3] *
                                        inv_three_metric[stage][c1].p->m[8] );
    //1,1 Component
    three_metric[stage][c1].p->m[4] = ( inv_three_metric[stage][c1].p->m[0] *
                                        inv_three_metric[stage][c1].p->m[8] -
                                        inv_three_metric[stage][c1].p->m[2] *
                                        inv_three_metric[stage][c1].p->m[6] );
    //1,2 Component
    three_metric[stage][c1].p->m[5] = ( inv_three_metric[stage][c1].p->m[2] *
                                        inv_three_metric[stage][c1].p->m[3] -
                                        inv_three_metric[stage][c1].p->m[0] *
                                        inv_three_metric[stage][c1].p->m[5] );
    //2,0 Component
    three_metric[stage][c1].p->m[6] = ( inv_three_metric[stage][c1].p->m[3] *
                                        inv_three_metric[stage][c1].p->m[7] -
                                        inv_three_metric[stage][c1].p->m[4] *
                                        inv_three_metric[stage][c1].p->m[6] );
    //2,1 Component
    three_metric[stage][c1].p->m[7] = ( inv_three_metric[stage][c1].p->m[1] *
                                        inv_three_metric[stage][c1].p->m[6] -
                                        inv_three_metric[stage][c1].p->m[0] *
                                        inv_three_metric[stage][c1].p->m[7] );
    //2,2 Component
    three_metric[stage][c1].p->m[8] = ( inv_three_metric[stage][c1].p->m[0] *
                                        inv_three_metric[stage][c1].p->m[4] -
                                        inv_three_metric[stage][c1].p->m[1] *
                                        inv_three_metric[stage][c1].p->m[3] );

    //finally divide by the determinant
    three_metric[stage][c1].ScalarMult(1.0/det);
   }
}


//******************************************************************************
//Name:  update_u_dns                                                          *
//                                                                             *
//Purpose:  returns the value of T_{\mu\nu} at a specified grid point          *
//                                                                             *
//Takes: an int telling what stage of the RK and 3 ints telling the grid       *
//       location                                                              *
//******************************************************************************
void Particle::update_u_dns(int stage)
{
  for(c1 = 0; c1 < numParticles; c1++)
   {
     //first execution of the program
     if( last_stage == -1 )
      {
        T_temp1x3 <= three_metric[stage][c1].
                      Contract((u_up[stage][c1] - shift[stage][c1]),1,0);

        T_temp1x1 <= (u_up[stage][c1] - shift[stage][c1]).
                       Contract(T_temp1x3,0,0);

        temp = 1.0 / sqrt(   lapse[stage][c1]*lapse[stage][c1]
                           - T_temp1x1.p->m[0] );

        u_dn[stage][c1] <= T_temp1x3 * temp;
     }

     T_temp1x3 <= inv_three_metric[stage][c1].Contract(u_dn[stage][c1],0,0);
     T_temp1x1 <= u_dn[stage][c1].Contract(T_temp1x3,0,0);

     u0_dn[stage][c1] = (shift[stage][c1].Contract(u_dn[stage][c1],0,0)).p->m[0]
                          - lapse[stage][c1]*sqrt(1 + T_temp1x1.p->m[0]);
   }
}

//******************************************************************************
//Name:  update_u_ups                                                          *
//                                                                             *
//Purpose:  returns the value of T_{\mu\nu} at a specified grid point          *
//                                                                             *
//Takes: an int telling what stage of the RK and 3 ints telling the grid       *
//       location                                                              *
//******************************************************************************
void Particle::update_u_ups(int stage)
{

  for(c1 = 0; c1 < numParticles; c1++)
   {
    //0 Component of contravariant four-velocity
    u0_up[stage][c1] =  inv_g[stage][c1].p->m[0] * u0_dn[stage][c1]
                      + inv_g[stage][c1].p->m[1] * u_dn[stage][c1].p->m[0]
                      + inv_g[stage][c1].p->m[2] * u_dn[stage][c1].p->m[1]
                      + inv_g[stage][c1].p->m[3] * u_dn[stage][c1].p->m[2];

    //0 Component of spatial portion of contravariant four-velocity
    u_up[stage][c1].p->m[0] =  inv_g[stage][c1].p->m[4] * u0_dn[stage][c1]
                             + inv_g[stage][c1].p->m[5] * u_dn[stage][c1].p->m[0]
                             + inv_g[stage][c1].p->m[6] * u_dn[stage][c1].p->m[1]
                             + inv_g[stage][c1].p->m[7] * u_dn[stage][c1].p->m[2];

    //1 Component of spatial portion of contravariant four-velocity
    u_up[stage][c1].p->m[1] =  inv_g[stage][c1].p->m[8]  * u0_dn[stage][c1]
                             + inv_g[stage][c1].p->m[9]  * u_dn[stage][c1].p->m[0]
                             + inv_g[stage][c1].p->m[10] * u_dn[stage][c1].p->m[1]
                             + inv_g[stage][c1].p->m[11] * u_dn[stage][c1].p->m[2];

    //2 Component of spatial portion of contravariant four-velocity
    u_up[stage][c1].p->m[2] =  inv_g[stage][c1].p->m[12] * u0_dn[stage][c1]
                             + inv_g[stage][c1].p->m[13] * u_dn[stage][c1].p->m[0]
                             + inv_g[stage][c1].p->m[14] * u_dn[stage][c1].p->m[1]
                             + inv_g[stage][c1].p->m[15] * u_dn[stage][c1].p->m[2];
    //Finally pack the 4-dimensional arrays
    four_U_dn[stage][c1].p->m[0] = u0_dn[stage][c1];
    four_U_dn[stage][c1].p->m[1] = u_dn[stage][c1].p->m[0];
    four_U_dn[stage][c1].p->m[2] = u_dn[stage][c1].p->m[1];
    four_U_dn[stage][c1].p->m[3] = u_dn[stage][c1].p->m[2];

    four_U_up[stage][c1].p->m[0] = u0_up[stage][c1];
    four_U_up[stage][c1].p->m[1] = u_up[stage][c1].p->m[0];
    four_U_up[stage][c1].p->m[2] = u_up[stage][c1].p->m[1];
    four_U_up[stage][c1].p->m[3] = u_up[stage][c1].p->m[2];

   }

}


//******************************************************************************
//Name:  RK_step                                                               *
//                                                                             *
//Purpose:                                                                     *
//                                                                             *
//Takes:                                                                       *
//                                                                             *
//******************************************************************************
void Particle::RK_step(stage)
{
 if( stage == 1 ) { next_stage = 0; time_step = delta_t / 2.0; }
 if( stage == 2 ) { next_stage = 1; time_step = 0.0; }
 if( stage == 0 ) { next_stage = 2; time_step = delta_t; }

 for(c1 = 0; c1 < numParticles; c1++)
  {
     temp = time_step / u0_up[stage][c1];
     pos[next_stage][c1] = pos[stage][c1] + u_up[stage][c1] * temp;

     //construct the contration (bilinear form) of the covariant 4-velocity
     //with the spatial derivative of the inverse metric
     //put the result in a temporary variable
     T_temp1x3 <= four_U_dn[stage][c1].Contract(
                  d_inv_g[stage][c1].Contract(four_U_dn[stage][c1],1,0),0,0);

     //form the RK step for the spatial covariant velocity
     u_dn[next_stage][c1] =   u_dn[stage][c1]
                            - 0.5 / u0_up[stage][c1] * time_step * T_temp1x3
                            + u0_dn[stage][c1] / norm[stage][c1]
                              * d_norm[stage][c1] * time_step;
  }

}

//******************************************************************************
//Name:  det_g                                                                 *
//                                                                             *
//Purpose:  calculates the determinant of the metric at a grid point           *
//                                                                             *
//Takes: an int telling what stage of the RK and 3 ints telling the grid       *
//       location                                                              *
//******************************************************************************
double  Particle::det_g(int stage, int i, int j, int k)
{
        //form the determinant of the Jacobian
//        detJ =   ( J[1][1] * J[2][2] - J[1][2] * J[2][1] ) * J[0][0]
//               + ( J[1][2] * J[2][0] - J[1][0] * J[2][2] ) * J[0][1]
//               + ( J[1][0] * J[2][1] - J[1][1] * J[2][0] ) * J[0][2];
//        detJ_inv = 1.0 / detJ;
 return -999.9999;
}

//******************************************************************************
//Name:  GetIndices                                                            *
//                                                                             *
//Purpose:  get the indices for the closest grid point                         *
//                                                                             *
//Takes: pointer to a position array and a pointer to an index array           *
//******************************************************************************
void Particle::GetIndices(double *position, int *indices)
{


//set the index array
	indices[0] = (int) ( (position[0] - spatial_extent[X_MIN])           /
		                  (spatial_extent[X_MAX] - spatial_extent[X_MIN])  )
                        * (number_of_points[0] - 1);

	indices[1] = (int) ( (position[1] - spatial_extent[Y_MIN])           /
		                  (spatial_extent[Y_MAX] - spatial_extent[Y_MIN])  )
                        * (number_of_points[1] - 1);

	indices[2] = (int) ( (position[2] - spatial_extent[Z_MIN])           /
		                  (spatial_extent[Z_MAX] - spatial_extent[Z_MIN])  )
                        * (number_of_points[2] - 1);
}

//******************************************************************************
//Name:  GetIndices                                                            *
//                                                                             *
//Purpose:  get the indices for the closest grid point                         *
//                                                                             *
//Takes: pointer to a position array and a pointer to a double index array     *
//******************************************************************************
void Particle::GetIndices(double *position, double *indices)
{

//set the index array
	indices[0] = (position[0] - spatial_extent[X_MIN])            /
		          (spatial_extent[X_MAX] - spatial_extent[X_MIN])
                * (number_of_points[0] - 1);

	indices[1] = (position[1] - spatial_extent[Y_MIN])            /
		          (spatial_extent[Y_MAX] - spatial_extent[Y_MIN])
                * (number_of_points[1] - 1);


	indices[2] = (position[2] - spatial_extent[Z_MIN])            /
		          (spatial_extent[Z_MAX] - spatial_extent[Z_MIN])
                * (number_of_points[2] - 1);
}

//******************************************************************************
//Name:  GetPosition                                                           *
//                                                                             *
//Purpose:  get the tensor containing the position of a grid based on its      *
//          indices                                                            *
//                                                                             *
//Takes: pointer to an index array                                             *
//******************************************************************************
tensor Particle::GetPosition(int *indices)
{
   tensor pos;

	int i;


//fill the pos tensor
	for(i =	0; i < 3; i++)
	{
		pos.Set((double)indices[i] * Deltas[i] + spatial_extent[2*i], i);
	}

   return pos;
}

//******************************************************************************
//Name:  GetPosition                                                           *
//                                                                             *
//Purpose:  get the position of a grid point based on it indices               *
//                                                                             *
//Takes: pointer to an index array and a pointer to a position array           *
//******************************************************************************
void Particle::GetPosition(int *indices, double *position)
{

	int i;


//fill the position array
	for(i =	0; i < 3; i++)
	{
		position[i] = (double)indices[i] * Deltas[i] + spatial_extent[2*i];
	}
}

//******************************************************************************
//Name:  SetSmoother                                                           *
//                                                                             *
//Purpose:  set the smoothing method for the particle or the grid              *
//                                                                             *
//Takes: int id and an int for the smoothing method                            *
//******************************************************************************
void Particle::SetSmoother(int id, int smoothing_method)
{
   if(id == PARTICLE_SMOOTHING)
	   particle_smoothing_method = smoothing_method;

   if(id == GRID_SMOOTHING)
	   grid_smoothing_method = smoothing_method;
}

//******************************************************************************
//Name:  SetGridRadius                                                         *
//                                                                             *
//Purpose:  set the smoothing radius for grid->particle                        *
//                                                                             *
//Takes: int that is the number of grid spacings                               *
//******************************************************************************
void Particle::SetGridRadius(int grid_rad)
{
   grid_radius = grid_rad;
   if( grid_radius == 0 )
    {
      grid_lower_bound = 0;
      grid_upper_bound = 2;
      smoothing_length = 1.0;
    }
   else
    {
      //the assymetry here is due to the for loop convention of
      //for(i = grid_lower_bound; i < grid_upper_bound; i++)
      //which could also be written as
      //for(i = grid_lower_bound; i < grid_upper_bound + 1; i++)
      //using a symmetric definition
      grid_lower_bound = -grid_radius - 1;
      grid_upper_bound = grid_radius + 2;
      smoothing_length = grid_radius;
    }
}



//******************************************************************************
//Name:  Dist                                                                  *
//                                                                             *
//Purpose:  get the Euclidean distance associated with a position              *
//                                                                             *
//Takes: pointer to the position array                                         *
//******************************************************************************
double Particle::Dist(double *position)
{
	double dist_ret;

    dist_ret = sqrt( position[0]*position[0] +
					      position[1]*position[1] +
					      position[2]*position[2] );

    return dist_ret;
}


//******************************************************************************
//Name:  GetInterpTensor                                                       *
//                                                                             *
//Purpose:  get the interpolated tensor at an arbitrary position               *
//                                                                             *
//Takes: pointer to the position as packed in a tensor                         *
//******************************************************************************
tensor Particle::GetInterpTensor(tensor &pos)
{
   double position[3];
   tensor ret_val;

   for(int i = 0; i < 3; i++)  position[i] = pos.Val(i);

   ret_val <= GetInterpTensor(position);

   return ret_val;
}

//******************************************************************************
//Name:  GetInterpTensor                                                       *
//                                                                             *
//Purpose:  get the interpolated tensor at an arbitrary position               *
//                                                                             *
//Takes: pointer to the position array                                         *
//******************************************************************************
tensor Particle::GetInterpTensor(double *position)
{
 int i, j, k, m, n;
 double scaled_pos[3];
 int grid_pos[3];
 double W;

 //initialize member data tensors
 sum_Wm   <= 0.0*sum_Wm;

//Set the reference position
  for(i = 0; i < 3; i++) reference_pos[i] = position[i];

//Get the real indices for the scaled position
  GetIndices(reference_pos, scaled_pos);

//Create the return tensor and set it to zero
// ret_dummy <= 0.0*GetTensor(reference_pos);

//Reset the normalization
  normalization = 0.0;

//step through the floor steping
  for( i = grid_lower_bound; i < grid_upper_bound; i++)
   {
	 for( j = grid_lower_bound; j < grid_upper_bound; j++)
	  {
	   for( k = grid_lower_bound; k < grid_upper_bound; k++)
		 {
         //Get the indices of the lowest corner of the cell
			grid_pos[0] = (int)floor( scaled_pos[0] + i );
			grid_pos[1] = (int)floor( scaled_pos[1] + j );
         grid_pos[2] = (int)floor( scaled_pos[2] + k );

         //find the position this grid point
         GetPosition(grid_pos,current_pos);

         //form the difference array and its magnitude (for use by smoother)
         //and the difference position
         y = 0;
         for( m = 0; m < 3; m++ )
			 {
           pos_diff[m] = reference_pos[m] - current_pos[m];
           y += pos_diff[m]*pos_diff[m];
           d.Set(pos_diff[m], m);
			 }
         y = sqrt(y)/smoothing_length;
         if( y > 1.0 + SMOOTHING_TOL ) continue;

         //get the (psuedo) scalar smoothing value
         W = Smoother(CWM_SMOOTHER);

         //get the derivative of the smoothing
         Smoother_Deriv(CWM_SMOOTHER);

         //keep running sum
         sum_Wm <= sum_Wm + Wm;
	    }
	  }
   }

//Normalize the interpolated tensor

 return ret_dummy;

}


//******************************************************************************
//Name:  Smoother                                                              *
//                                                                             *
//Purpose:  general function that mediates the various smoothings              *
//                                                                             *
//Takes: an int telling what type of smoothing                                 *
//       (grid->particle or particle->grid)                                    *
//******************************************************************************
double Particle::Smoother(int type)
{
   double W;

   switch(type)
   {
     case CWM_SMOOTHER:
                if( y > 1.0 + SMOOTHING_TOL )
                 {
                   W = 0.0;
                 }
                else
                 {
                   W = pow((1.0 - y*y), 3);
                 }
                break;

     case DEFAULT_SMOOTHING:
                W = pow((1.0 - y*y), 3);
                break;
   }
   normalization += fabs(W);
	return W;
}

//******************************************************************************
//Name:  Smoother_Deriv                                                        *
//                                                                             *
//Purpose:  general function that models the derivatives of the various        *
//          smoothings in Smoother                                            *
//                                                                             *
//Takes: an int telling what type of smoothing                                 *
//       (grid->particle or particle->grid)                                    *
//******************************************************************************
void Particle::Smoother_Deriv(int type)
{
   double q;
   int i;
   double arg;


  switch(type)
   {
     case CWM_SMOOTHER:
                if( y > 1.0 + SMOOTHING_TOL )
                 {
                   Wm <= 0.0*Wm;
                 }
                else
                 {
                   arg = -6.0*(1.0 - y*y)*(1.0 - y*y)/smoothing_length;
                   Wm.p->m[0] = arg*pos_diff[0];
                   Wm.p->m[1] = arg*pos_diff[1];
                   Wm.p->m[2] = arg*pos_diff[2];
                 }
                break;

     case DEFAULT_SMOOTHING:
                break;
   }
}


//******************************************************************************
//Name:  SetParams                                                             *
//                                                                             *
//Purpose:  to make the simulation parameters known to the particle class      *
//                                                                             *
//Takes: a Params structure                                                    *
//******************************************************************************
void Particle::SetParams(Param simParamsIn)
{
	masses[0] = simParamsIn.mass0;
	masses[1] = simParamsIn.mass1;

	pos[1][0].p->m[0]  = simParamsIn.particleCenter0.p->m[1];
	pos[1][0].p->m[1]  = simParamsIn.particleCenter0.p->m[2];
	pos[1][0].p->m[2]  = simParamsIn.particleCenter0.p->m[3];
   u_up[1][0].p->m[0] = simParamsIn.particleVelocity0.p->m[1];
   u_up[1][0].p->m[1] = simParamsIn.particleVelocity0.p->m[2];
   u_up[1][0].p->m[2] = simParamsIn.particleVelocity0.p->m[3];
	pos[1][1].p->m[0]  = simParamsIn.particleCenter1.p->m[1];
	pos[1][1].p->m[1]  = simParamsIn.particleCenter1.p->m[2];
	pos[1][1].p->m[2]  = simParamsIn.particleCenter1.p->m[3];
   u_up[1][1].p->m[0] = simParamsIn.particleVelocity1.p->m[1];
   u_up[1][1].p->m[1] = simParamsIn.particleVelocity1.p->m[2];
   u_up[1][1].p->m[2] = simParamsIn.particleVelocity1.p->m[3];
   spatial_extent[X_MIN] = simParamsIn.minimums.Val(0);
   spatial_extent[Y_MIN] = simParamsIn.minimums.Val(1);
   spatial_extent[Z_MIN] = simParamsIn.minimums.Val(2);
   spatial_extent[X_MAX] = simParamsIn.maximums.Val(0);
   spatial_extent[Y_MAX] = simParamsIn.maximums.Val(1);
   spatial_extent[Z_MAX] = simParamsIn.maximums.Val(2);
   Deltas[X_ID] = ( spatial_extent[X_MAX] - spatial_extent[X_MIN] )
                  / ( simParamsIn.max.Val(X_ID + 1) - 1 );
   Deltas[Y_ID] = ( spatial_extent[Y_MAX] - spatial_extent[Y_MIN] )
                  / ( simParamsIn.max.Val(Y_ID + 1) - 1 );
   Deltas[Z_ID] = ( spatial_extent[Z_MAX] - spatial_extent[Z_MIN] )
                  / ( simParamsIn.max.Val(Z_ID + 1) - 1 );
   Delta_mag = sqrt(    Deltas[X_ID]*Deltas[X_ID]
                     +  Deltas[Y_ID]*Deltas[Y_ID]
                     +  Deltas[Z_ID]*Deltas[Z_ID] );
   inv_Delta_mag = 1.0 / Delta_mag;
   inv_Deltas[X_ID] = 1.0 / Deltas[X_ID];
   inv_Deltas[Y_ID] = 1.0 / Deltas[Y_ID];
   inv_Deltas[Z_ID] = 1.0 / Deltas[Z_ID];
   number_of_points[0] = simParamsIn.max.Val(1);
   number_of_points[1] = simParamsIn.max.Val(2);
   number_of_points[2] = simParamsIn.max.Val(3);
//   smoothing_length    = simParamsIn.smoothLength0.Val(1);
   SetGridRadius(0);
   delta_t = simParamsIn.delta.p->m[0];

}




