void main(void);
#include "C:\My Documents\Programming\tensor++\tensor.hpp"

void main(void)
{
   tensor A(3,2,2,2), B(2,2,2), C(2,2,2), Q;
   int i, j, k;
   char *dump;
   FILE* dude;

   dude = fopen("contract_test.txt","w");
   dump = NULL;

   for(i = 0; i < 2; i++)
     for(j = 0; j < 2; j++)
       for(k = 0; k < 2; k++)
         A.Set((i+1.0) - 2.0*(j+1.0) + (k+1.0),i,j,k);

  A.print("A",&dump);
  fprintf(dude,"A = \n");
  fprintf(dude,"%s",dump);

  for(i = 0; i < 2; i++)
    for(j = 0; j < 2; j++)
    {
       B.Set(     (i+1.0) + 2.0*(j+1.0),i,j);
       C.Set( 3.0*(j+1.0) -     (i+1.0),i,j);
    }

  B.print("B",&dump);
  fprintf(dude,"B=\n");
  fprintf(dude,"%s",dump);

  C.print("C",&dump);
  fprintf(dude,"C=\n");
  fprintf(dude,"%s",dump);

  Q <= (A.Contract(C,1,0)).Contract(B,0,1);
  Q.print("Q",&dump);
  fprintf(dude,"Q <= (A.Contract(C,1,0)).Contract(B,0,1)\n");
  fprintf(dude,"%s",dump);

  Q <= (A.Contract(B,0,1)).Contract(C,0,0);
  Q.print("Q",&dump);
  fprintf(dude,"Q <= (A.Contract(B,0,1)).Contract(C,0,0)\n");
  fprintf(dude,"%s",dump);

  Q <= (B.Contract(A,1,0)).Contract(C,1,0);
  Q.print("Q",&dump);
  fprintf(dude,"Q <= (B.Contract(A,1,0)).Contract(C,1,0)\n");
  fprintf(dude,"%s",dump);

  Q <= B.Contract(C.Contract(A,0,1),1,1);
  Q.print("Q",&dump);
  fprintf(dude,"Q <= B.Contract(C.Contract(A,0,1),1,1)\n");
  fprintf(dude,"%s",dump);

  Q <= (C.Contract(A,0,1)).Contract(B,1,1);
  Q.print("Q",&dump);
  fprintf(dude,"Q <= (C.Contract(A,0,1)).Contract(B,0,1)\n");
  fprintf(dude,"%s",dump);

  Q <= C.Contract(B.Contract(A,1,0),0,1);
  Q.print("Q",&dump);
  fprintf(dude,"Q <= C.Contract(B.Contract(A,1,0),0,1)\n");
  fprintf(dude,"%s",dump);

  fclose(dude);
}
