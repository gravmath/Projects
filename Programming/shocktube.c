/* $Header$ */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>



#define NPARTICLES 220
#define TRIAG_KERNEL 1

/* smoothing length */
#define HSML 0.095

/* artificial viscosity coefficient */
const double beta = 1.0;
const double epsilon = 0.01;

/* number of time steps */
const int nsteps = 200;

/* time step size */
const double dt = 0.005;

/* pi */
#define PI_CONST 3.14159265358979324

struct particle {
  double x, v, m, vol, eps;
  double xdot, vdot, voldot, epsdot;
  double p;
  int fixed;
};



int step;
double time;

struct particle particles[NPARTICLES];


#if TRIAG_KERNEL
double kernel (const double dx)
{
  if (fabs(dx) >= HSML) {
    return 0;
  } else {
    return (HSML - fabs(dx)) / (HSML*HSML);
  }
}

double grad_kernel (const double dx)
{
  if (fabs(dx) >= HSML) {
    return 0;
  } else if (fabs(dx) < 1e-10) {
    return 0;
  } else {
    return - (dx<0 ? -1 : +1) / (HSML*HSML);
  }
}

#else
double kernel (const double dx)
{
  static double coeff = 1.0/(2.0*HSML);
  if (fabs(dx) > 3.0*HSML) {
    return 0.0;
  } else {
    return coeff * exp(-fabs(dx)/HSML);
  }
}

double grad_kernel (const double dx)
{
  static double coeff = -1.0/(2.0*HSML*HSML);
  if (fabs(dx) > 3.0*HSML) {
    return 0.0;
  } else {
    return (dx < 0 ? -1 : 1) * coeff * exp(-fabs(dx)/HSML);
  }
}
#endif

double potential (const double x, const double rho)
{
  /* Newtonian gravitational potential density for central mass = 100 */
  return (100.0 * rho)/x;
}

double grad_potential (const double x, const double rho) 
{
  /* d/dx Newtonian gravitational potential density for central mass = 100 */
  return -(100.0 * rho)/(x*x);
}

void output (void)
{
  static int firsttime = 1;
  FILE * file;
  int i;
  file = fopen ("shocktube.out", firsttime ? "w" : "a");
  if (firsttime) {
    fprintf (file, "# n t   i   x v m vol eps   fixed   p   xdot vdol voldot epsdot\n");
  }
  for (i=0; i<NPARTICLES; ++i) {
    fprintf (file, "%d %g   %d   %g %g %g %g %g   %d   %g   %g %g %g %g\n",
	     step,
	     time,
	     i,
	     particles[i].x,
	     particles[i].v,
	     particles[i].m,
	     particles[i].vol,
	     particles[i].eps,
	     particles[i].fixed,
	     particles[i].p,
	     particles[i].xdot,
	     particles[i].vdot,
	     particles[i].voldot,
	     particles[i].epsdot);
  }
//  fprintf (file, "\n");
  fclose (file);
  firsttime = 0;
}



void initialise (void)
{
  int i;
  step = 0;
  time = 0.0;
  for (i=0; i<NPARTICLES; ++i) {
    particles[i].x = (i - NPARTICLES/2 + 0.5) * 0.01;
    particles[i].v = 0.0;
    particles[i].m = i<NPARTICLES/2 ? 0.02 : 0.01;
    particles[i].vol = 0.01;
    particles[i].eps = 1.0;
    particles[i].fixed = i<10 || i>=210;
  }
}



void statistics (void)
{
  double V, M, P, E, CM;
  int i;
  V = M = P = E = CM = 0;
  for (i=0; i<NPARTICLES; ++i) {
    if (! particles[i].fixed) {
      V += particles[i].vol;
      M += particles[i].m;
      P += particles[i].m * particles[i].v;
      E += (0.5 * particles[i].m * particles[i].v * particles[i].v
	    + particles[i].m * particles[i].eps);
      CM += particles[i].m * particles[i].x;
    }
  }
  CM /= M;
  printf ("current time:   %g\n", time);
  printf ("total mass:     %g\n", M);
  printf ("total volume:   %g\n", V);
  printf ("total momentum: %g\n", P);
  printf ("total energy:   %g\n", E);
  printf ("centre of mass: %g\n", CM);
  printf ("\n");
}


/* pressure of particle i */
double pressure (const int i)
{
  return particles[i].m / particles[i].vol * particles[i].eps;
}

/* time derivatives of particle i that are caused by particle j */
double vdot (const int i, const int j)
{
  /* - 1/rho grad p */
  return (- 1/particles[i].m
	  * particles[i].vol * particles[j].vol
	  * (particles[j].p + particles[i].p)
	  * grad_kernel(particles[i].x - particles[j].x));
}

double voldot (const int i, const int j)
{
  /* V div v */
  return (particles[i].vol * particles[j].vol
	  * (particles[j].v - particles[i].v)
	  * grad_kernel(particles[i].x - particles[j].x));
}

double epsdot (const int i, const int j)
{
  /* - p/rho div v */
  return (- particles[i].p / particles[i].m
	  * particles[i].vol * particles[j].vol
	  * (particles[j].v - particles[i].v)
	  * grad_kernel(particles[i].x - particles[j].x));
}

double artvisc_vdot (const int i, const int j)
{
  const double dx = particles[i].x - particles[j].x;
  const double dv = particles[i].v - particles[j].v;
  const double rhoij = 0.5 * (particles[i].m / particles[i].vol
			      + particles[j].m / particles[j].vol);
  const double dq = dx / HSML;
  const double tmp = (dq * dv) / (dq*dq + epsilon);
  const double pi = beta * tmp * tmp / rhoij;
  if (tmp >= 0) {
    return 0;
  } else {
    return (- particles[j].m * pi
	    * grad_kernel(particles[i].x - particles[j].x));
  }
}

double artvisc_epsdot (const int i, const int j)
{
  const double dx = particles[i].x - particles[j].x;
  const double dv = particles[i].v - particles[j].v;
  const double rhoij = 0.5 * (particles[i].m / particles[i].vol
			      + particles[j].m / particles[j].vol);
  const double dq = dx / HSML;
  const double tmp = (dq * dv) / (dq*dq + epsilon);
  const double pi = beta * tmp * tmp / rhoij;
  if (tmp >= 0) {
    return 0;
  } else {
    return (- particles[j].m
	    * 0.5 * pi
	    * (particles[j].v - particles[i].v)
	    * grad_kernel(particles[i].x - particles[j].x));
  }
}



void rhs (void)
{
  int i, j;
  
  for (i=0; i<NPARTICLES; ++i) {
    particles[i].p = pressure(i);
  }
  
  for (i=0; i<NPARTICLES; ++i) {
    particles[i].xdot = particles[i].v;
    
    particles[i].vdot = 0;
    particles[i].voldot = 0;
    particles[i].epsdot = 0;
    for (j=0; j<NPARTICLES; ++j) {
      particles[i].vdot += vdot(i,j) + artvisc_vdot(i,j);
      particles[i].voldot += voldot(i,j);
      particles[i].epsdot += epsdot(i,j) + artvisc_epsdot(i,j);
    }
  }
}



void euler (void)
{
  int i;
  
  rhs ();
  
  for (i=0; i<NPARTICLES; ++i) {
    if (! particles[i].fixed) {
      particles[i].x += dt * particles[i].xdot;
      particles[i].v += dt * particles[i].vdot;
      particles[i].vol += dt * particles[i].voldot;
      particles[i].eps += dt * particles[i].epsdot;
    }
  }
  
  ++step;
  time += dt;
}



int main (const int argc, char ** const argv)
{
  int n;
  printf ("shocktube\n");
  printf ("initial data\n");
  initialise ();
  rhs ();
  statistics ();
  output ();
  for (n=0; n<nsteps; ++n) {
    printf ("time step %d\n", n+1);
    euler ();
    rhs ();
    statistics ();
    output ();
  }
  printf ("done.\n");
  return 0;
}
