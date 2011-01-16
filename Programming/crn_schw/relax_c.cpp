#include <stdio.h>
#include <math.h>
#include "array.hpp"
#include <dos.h>

struct time _time1;
struct time _time2;


void main(void);
#define LATTICE_SIZE  35
#define OUTPUT_FREQ   50000
#define LOOP          50000

void main(void)
{
   gettime(&_time1);

	//Define Variables
   array lattice(3,LATTICE_SIZE,LATTICE_SIZE, LATTICE_SIZE);
	int i,j, k, counter;
   double coeff_1_6 = 1.0/6.0;

	//perform desired iterations
	for( counter = 0; counter < LOOP; counter++)
	{
	  //relax the lattice
	  for( i = 1; i < LATTICE_SIZE - 1; i++)
		for( j = 1; j < LATTICE_SIZE - 1; j++)
        for( k = 1; k < LATTICE_SIZE - 1; k++)
        	{
		   	lattice.Set( coeff_1_6*(
                                     lattice.Val(i+1,j,k) +
                                     lattice.Val(i-1,j,k) +
                                     lattice.Val(i,j+1,k) +
                                     lattice.Val(i,j-1,k) +
                                     lattice.Val(i,j,k-1) +
                                     lattice.Val(i,j,k+1) ), i,j,k);
		   }

	}

  gettime(&_time2);

  int temp =  (_time2.ti_hour - _time1.ti_hour)*3600
            + (_time2.ti_min - _time1.ti_min)*60
            + (_time2.ti_sec - _time1.ti_sec);
  printf("The time was %d secs",temp);

}