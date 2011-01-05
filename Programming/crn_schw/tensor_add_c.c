#include <stdio.h>
#include <math.h>
#include <dos.h>

struct time _time1;
struct time _time2;

#define SIZE 30

void main(void);

void main(void)
{
   gettime(&_time1);
   double A[SIZE][SIZE][SIZE], B[SIZE][SIZE][SIZE];
   double tensor3[SIZE][SIZE][SIZE];
   int i, j, k, counter;

   for(i = 0; i < SIZE; i++)
    for(j = 0; j < SIZE; j++)
     for(k = 0; k < SIZE; k++)
      {
       A[i][j][k] = 0;
       B[i][j][k] = 0;
       tensor3[i][j][k] = 0;
      }

   for( counter = 0; counter < 50001; counter ++)
    {
       for(i = 0; i < SIZE; i++)
        for(j = 0; j < SIZE; j++)
         for(k = 0; k < SIZE; k++)
          {
           tensor3[i][j][k] = A[i][j][k] + 3*B[i][j][k];
          }
    }
  gettime(&_time2);

  int temp = (_time2.ti_min - _time1.ti_min)*60 + (_time2.ti_sec - _time1.ti_sec);
  printf("The time was %d secs",temp);


}

