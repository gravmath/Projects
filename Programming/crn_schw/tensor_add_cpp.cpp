#include <stdio.h>
#include <math.h>
#include <dos.h>

#include "tensor.hpp"

struct time _time1;
struct time _time2;

#define SIZE 30
#define LOOP 5000

void main(void);

void main(void)
{
   gettime(&_time1);
   tensor A(3,SIZE,SIZE,SIZE), B(3,SIZE,SIZE,SIZE);
   tensor tensor3(3,SIZE,SIZE,SIZE);
//   int i, j, k;
   int counter;

   for( counter = 0; counter < LOOP; counter ++)
    {
      tensor3 = A + 3*B;
    }
  gettime(&_time2);

  int temp = (_time2.ti_min - _time1.ti_min)*60 + (_time2.ti_sec - _time1.ti_sec);
  printf("The time was %d secs",temp);


}

