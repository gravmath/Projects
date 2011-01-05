/*This is a test of the array class component of the tensor++
formalism.  The test is as follows:

Name:	test1

Date: 4/5/03

Purpose:  Test the basic member function of array

*/

#include "C:\My Documents\Programming\tensor++\array.hpp"
void main(void);

void main(void)
{
  array A(1,10), B(2,2,2);
  array C, D;
  char *dump;
  dump = NULL;
  FILE *skippy;
  skippy = fopen("array_test1.txt","w");

  for( int i = 0; i < 10; i++)
    A.Set(i*i,i);

  A.SetName("A");
  A.print("A",&dump);
  fprintf(skippy,"A.Set(i*i,i)\n");
  fprintf(skippy,"%s",dump);
  fprintf(skippy,"A.max() = %g\n\n",A.max());

  B.SetName("B");
  B.Set(-99,0,0);
  B.Set(-99,0,1);
  B.Set(-99,1,0);
  B.Set(-99,1,1);
  B.print("B",&dump);
  fprintf(skippy,"Contents of B\n");
  fprintf(skippy,"%s\n\n",dump);

  C <= B;
  D <= A;
  C.print("C",&dump);
  fprintf(skippy,"Contents of C and D\n");
  fprintf(skippy,"%s",dump);
  D.print("D",&dump);
  fprintf(skippy,"%s",dump);

  C.Resize(1,5);
  C.print("C-resized",&dump);
  fprintf(skippy,"Contents of C after Resize\n");
  fprintf(skippy,"%s",dump);
}
