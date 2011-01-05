#include <stdio.h>
#include <math.h>
void main(void);

void main(void)
{
 	double **my_array;

   my_array = new double* [3];

   for( int i = 0; i < 3; i++ )
   {
      my_array[i] = new double[2];
   }

   for( int i = 0; i < 3; i++)
     for( int j = 0; j < 2; j++)
       my_array[i][j] = 2*i+j;


}
