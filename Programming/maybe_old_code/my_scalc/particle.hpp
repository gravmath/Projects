//******************************************************************************
//Name:  Particle.hpp                                                          *
//                                                                             *
//Purpose:  header file for defining the particle class                        *
//                                                                             *
//Modification History:  12/28/98 - first draft                                *
// *****************************************************************************
#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "tensor.hpp"
#include "param.hpp"
#include "MSField.hpp"

#define X_MIN 0
#define X_MAX 1
#define Y_MIN 2
#define Y_MAX 3
#define Z_MIN 4
#define Z_MAX 5

#define X_ID 0
#define Y_ID 1
#define Z_ID 2

#define SPATIAL_TOL  1e-10

#define PARTICLE_SMOOTHING 0
#define GRID_SMOOTHING     1

#define CWM_SMOOTHER 0
#define DEFAULT_SMOOTHING 9999
#define SMOOTHING_TOL 1e-7
class MSField;
class Particle {

public:
	//Constructor
	Particle();
	Particle(int number_of_particles);

	//Destructor
	~Particle();


	// Set the various parameters for a particle
	void 		SetPos					(int    particle_label,
                                  int    stage,
			  		            	    double x,
         					   		 double y,
             	 	 	 				 double z);
	void 		SetVel					(int    particle_label,
                                  int    stage,
            		   			 	 double vx,
			     			           	 double vy,
         					       	 double vz);
	void 		SetMass					(int    particle_label,
		         						 double mass);
	void 		SetSmoothingLength	(int    particle_label,
		            		          double length);
   void 		SetTime					(double time);


	// Return the various parameters for a particle
	tensor 	GetPos					(int particle_label, int stage);
	tensor 	GetVel					(int particle_label, int stage);
	double 	GetPos					(int particle_label,
                                  int stage,
		           						 int component);
	double 	GetVel					(int particle_label,
                                  int stage,
		           						 int component);
	double 	GetMass					(int particle_label);
	double 	GetSmoothingLength	(int particle_label);
   double 	GetTime					(void);

   void    	SetSpatialExtent		(int id, double value);
   void    	SetSpatialExtent		(double *value);
   void    	GetIndices				(double *position, int *indices);
   void    	GetIndices				(double *position, double *indices);
   void    	GetPosition				(int *indices, double *position);
   tensor  	GetPosition				(int *indices);
   void    	SetSmoother				(int id, int smoothing_method);
   void    	SetGridRadius			(int grid_rad);
   tensor  	GetInterpTensor		(tensor &pos);
   tensor  	GetInterpTensor		(double *position);
   int     	GetNumIndices			(int id);

	//Hydrodynamics
   tensor 	matterVar				(int stage, int i, int j, int k);

   //Field
   void 		RegisterField			(MSField *field);

   //Simulation Parameters
   void 		SetParams				(Param simParamsIn);

private:
	//array of particles
   int    	numParticles;
	tensor 	*pos[3];                 //the '3' index refers to the RK stage
	tensor 	*u_dn[3];                //the '3' index refers to the RK stage
	tensor 	*inv_g[3];               //the '3' index refers to the RK stage
   tensor 	*u_up[3];                //the '3' index refers to the RK stage
   double   *u0_dn[3];					 //the '3' index refers to the RK stage
   double   *u0_up[3];               //the '3' index refers to the RK stage
   tensor   *d_inv_g[3];             //the '3' index refers to the RK stage
   tensor   *shift[3];               //the '3' index refers to the RK stage
   double   *lapse[3];               //the '3' index refers to the RK stage
   tensor   *three_metric[3];        //the '3' index refers to the RK stage
   tensor   *inv_three_metric[3];    //the '3' index refers to the RK stage
   tensor   *four_U_up[3];           //the '3' index refers to the RK stage
   tensor   *four_U_dn[3];           //the '3' index refers to the RK stage
   double   *norm[3];                //the '3' index refers to the RK stage
   tensor   *d_norm[3];              //the '3' index refers to the RK stage
   double 	inv_lapse_val;
   double   inv_lapse_sqr;
	double 	*masses;
   double 	*smoothing_lengths;
   char   	name[10];
   double   current_time;
   double 	last_update_time;
   int      next_stage;
   int      current_stage;
   int      last_stage;
   double   time_step;
   double   delta_t;
   double   det;

	//utility variables
	int 		c1;
   tensor   T_temp4x4x3;
	tensor 	T_temp4x4;
   tensor   T_temp3x3;
   tensor   T_temp1x3;
   tensor   T_temp1x1;
   tensor   inv_g_temp;
   double   temp;


   //field variables
   MSField *my_field_partner;

   //parameters related to the smoothing performed on the interpolated tensors
   //associated with a position
	double 	W;
	double 	y;
	double 	normalization;

   int     particle_smoothing_method;
   int     grid_smoothing_method;
   int     grid_radius;
   double  Smoother(int type);
   void    Smoother_Deriv(int type);
   double  reference_pos[3];
   double  current_pos[3];
   double  pos_diff[3];
   int     number_of_points[3];
   int     indices[3];

	//tensors defined in for use in GetInterpTensor herein declared in
	//the interest of speed
    tensor  Wm, sum_Wm;
    //parameters related to the time of the current spatial field slice

    int	   is_time_set;

    //parameters related to the spatial grid structure of the field
    double  spatial_extent[6];
    double  Deltas[3];
    double  inv_Deltas[3];
    double  Delta_mag;
    double  inv_Delta_mag;
    double  Dist(double *position);

    //parameters related to the smoothing performed on the interpolated tensors
    //associated with a position
    int     out_of_range;
    int     is_periodic;
    int     grid_lower_bound;
    int     grid_upper_bound;
    double  smoothing_length;
    double  factor;
    int     is_series_set;
    int     is_deriv_field;
    int     is_dep_field;
	//tensors defined in for use in GetInterpTensor herein declared in
	//the interest of speed
    tensor  *tmp;
    tensor  a2, sum_a2, sum_da1, d_sum_a2;
    tensor 	d, ret_dummy, temp0, temp1, temp2;
    tensor 	grid_pos, pos_diff_tensor;

    void    update_inv_g(int stage);
    void    update_u_dns(int stage);
    void    update_u_ups(int stage);
    void    ADM_decomposition(int stage);
    void    RK_step(stage);
    double  det_g(int stage, int i, int j, int k);
    Param   particle_parms;
};




//******************************************************************************
//     INLINES BEGIN HERE                                                      *
//******************************************************************************

// The destructor needs to return the memory used for the
//	dynamically-allocated arrays
inline Particle::~Particle()
{

}

//******************************************************************************
//Name:  particle SetPos                                                       *
//******************************************************************************
inline void Particle::SetPos(int particle_label, int stage,
                             double x, double y, double z)
{
 	pos[particle_label][stage].p->m[0] = x;
 	pos[particle_label][stage].p->m[1] = y;
  	pos[particle_label][stage].p->m[2] = z;
}

//******************************************************************************
//Name:  particle SetVel                                                       *
//******************************************************************************
inline void Particle::SetVel(int particle_label, int stage,
                             double vx,double vy, double vz)
{
 	u_dn[particle_label][stage].p->m[0] = vx;
 	u_dn[particle_label][stage].p->m[1] = vy;
  	u_dn[particle_label][stage].p->m[2] = vz;
}

//******************************************************************************
//Name:  particle SetMass                                                      *
//******************************************************************************
inline void Particle::SetMass(int particle_label, double mass)
{
  masses[particle_label] = mass;
}

//******************************************************************************
//Name:  particle SetSmoothingLength                                           *
//******************************************************************************
inline void Particle::SetSmoothingLength(int particle_label, double length)
{
   smoothing_lengths[particle_label] = length;
}

//******************************************************************************
//Name:  particle SetTime                                                      *
//******************************************************************************
inline void Particle::SetTime(double time)
{
    current_time = time;
}

//******************************************************************************
//Name:  particle GetPos                                                       *
//******************************************************************************
inline tensor Particle::GetPos(int particle_label, int stage)
{
  tensor temp;
  temp.p->m[0] = pos[particle_label][stage].p->m[0];
  temp.p->m[1] = pos[particle_label][stage].p->m[1];
  temp.p->m[2] = pos[particle_label][stage].p->m[2];

  return temp;
}

//******************************************************************************
//Name:  particle GetVel                                                       *
//******************************************************************************
inline tensor Particle::GetVel(int particle_label, int stage)
{
  tensor temp;
  temp.p->m[0] = u_dn[particle_label][stage].p->m[0];
  temp.p->m[1] = u_dn[particle_label][stage].p->m[1];
  temp.p->m[2] = u_dn[particle_label][stage].p->m[2];

  return temp;
}

//******************************************************************************
//Name:  particle GetPos                                                       *
//******************************************************************************
inline double Particle::GetPos(int particle_label, int stage, int component)
{
  return pos[particle_label][stage].p->m[component];
}

//******************************************************************************
//Name:  particle GetVel                                                       *
//******************************************************************************
inline double Particle::GetVel(int particle_label, int stage, int component)
{
  return u_dn[particle_label][stage].p->m[component];
}

//******************************************************************************
//Name:  particle GetMass                                                      *
//******************************************************************************
inline double Particle::GetMass(int particle_label)
{
  return masses[particle_label];
}

//******************************************************************************
//Name:  particle GetSmoothingLength                                           *
//******************************************************************************
inline double Particle::GetSmoothingLength(int particle_label)
{
  return smoothing_lengths[particle_label];
}

//******************************************************************************
//Name:  particle GetTime                                                      *
//******************************************************************************
inline double Particle::GetTime(void)
{
  return current_time;
}

#endif
