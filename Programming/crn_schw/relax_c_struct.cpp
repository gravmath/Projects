#include <stdio.h>
#include <math.h>
#include "array.hpp"
#include <dos.h>

struct time _time1;
struct time _time2;

struct c_array { int num_ind; double *m; int *ranges; int *scales;};

void main(void);
#define LATTICE_SIZE  35
#define OUTPUT_FREQ   50000
#define LOOP          50000

void main(void)
{
   gettime(&_time1);

	//Define Variables
   struct c_array lattice;
   lattice.num_ind = 3;
   lattice.m = (double *)calloc(LATTICE_SIZE*LATTICE_SIZE*LATTICE_SIZE,sizeof(double));
   lattice.ranges = (int *)calloc(lattice.num_ind, sizeof(int));
   lattice.scales = (int *)calloc(lattice.num_ind, sizeof(int));

   lattice.ranges[0] = LATTICE_SIZE;
   lattice.ranges[1] = LATTICE_SIZE;
   lattice.ranges[2] = LATTICE_SIZE;

   lattice.scales[0] = LATTICE_SIZE*LATTICE_SIZE;
   lattice.scales[1] = LATTICE_SIZE;
   lattice.scales[2] = 1;

	int i,j,k, counter;
   double coeff_1_6 = 1.0/6.0;

	//initialize the lattice
	for( i = 0; i < LATTICE_SIZE; i++)
		for( j = 0; j < LATTICE_SIZE; j++)
        for( k = 0;  k < LATTICE_SIZE; k++)
   			lattice.m[i*lattice.scales[0] + j*lattice.scales[1] + k*lattice.scales[2]] = 0;

	//perform LOOP iterations
	for( counter = 0; counter < LOOP; counter++)
	{
	  //relax the lattice
	  for( i = 1; i < LATTICE_SIZE - 1; i++)
		for( j = 1; j < LATTICE_SIZE - 1; j++)
       for( k = 1; k < LATTICE_SIZE - 1; k++)
		{
			lattice.m[i*lattice.scales[0] + j*lattice.scales[1] + k*lattice.scales[2]]
         = coeff_1_6*
         (
         lattice.m[(i+1)*lattice.scales[0] + j*lattice.scales[1] + k*lattice.scales[2]] +
         lattice.m[(i-1)*lattice.scales[0] + j*lattice.scales[1] + k*lattice.scales[2]] +
         lattice.m[i*lattice.scales[0] + (j+1)*lattice.scales[1] + k*lattice.scales[2]] +
         lattice.m[i*lattice.scales[0] + (j-1)*lattice.scales[1] + k*lattice.scales[2]] +
         lattice.m[i*lattice.scales[0] + j*lattice.scales[1] + (k+1)*lattice.scales[2]] +
         lattice.m[i*lattice.scales[0] + j*lattice.scales[1] + (k-1)*lattice.scales[2]]
         );
		}

	}

  gettime(&_time2);

  int temp =  (_time2.ti_hour - _time1.ti_hour)*3600
            + (_time2.ti_min - _time1.ti_min)*60
            + (_time2.ti_sec - _time1.ti_sec);
  printf("The time was %d secs",temp);

}