#include <stdio.h>
#include <math.h>
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
   double lattice[LATTICE_SIZE][LATTICE_SIZE][LATTICE_SIZE];
	int i,j,k, counter;
   double coeff_1_6 = 1.0/6.0;

	//initialize the lattice
	for( i = 0; i < LATTICE_SIZE; i++)
		for( j = 0; j < LATTICE_SIZE; j++)
        for( k = 0;  k < LATTICE_SIZE; k++)
   			lattice[i][j][k] = 0;

	//perform LOOP iterations
	for( counter = 0; counter < LOOP; counter++)
	{
	  //relax the lattice
	  for( i = 1; i < LATTICE_SIZE - 1; i++)
		for( j = 1; j < LATTICE_SIZE - 1; j++)
       for( k = 1; k < LATTICE_SIZE - 1; k++)
		{
			lattice[i][j][k] = coeff_1_6*( lattice[i+1][j][k] +
                                        lattice[i-1][j][k] +
                                        lattice[i][j+1][k] +
                                        lattice[i][j-1][k] +
                                        lattice[i][j][k-1] +
                                        lattice[i][j][k+1]);
		}

	}

  gettime(&_time2);

  int temp =  (_time2.ti_hour - _time1.ti_hour)*3600
            + (_time2.ti_min - _time1.ti_min)*60
            + (_time2.ti_sec - _time1.ti_sec);
  printf("The time was %d secs",temp);

}
