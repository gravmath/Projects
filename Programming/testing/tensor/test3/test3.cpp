#include "C:\My Documents\Programming\tensor++\tensor.hpp"

#define SIZE 72

void main(void);

void main(void)
{
   double***** C;
   int i, j, k, l, m;
   tensor*** D;


   C = new double**** [SIZE];

   for(i = 0; i < SIZE; i++)
     C[i] = new double*** [SIZE];


   for(i = 0; i < SIZE; i++)
     for(j = 0; j < SIZE; j++)
       C[i][j] = new double** [SIZE];


   for(i = 0; i < SIZE; i++)
     for(j = 0; j < SIZE; j++)
       for(k = 0; k < SIZE; k++)
         C[i][j][k] = new double* [2];

   for(i = 0; i < SIZE; i++)
     for(j = 0; j < SIZE; j++)
       for(k = 0; k < SIZE; k++)
         for(l = 0; l < 2; l++)
         C[i][j][k][l] = new double [2];

   for(i = 0; i < SIZE; i++)
     for(j = 0; j < SIZE; j++)
       for(k = 0; k < SIZE; k++)
         for(l = 0; l < 2; l++)
           for(m = 0; m < 2; m++)
             C[i][j][k][l][m] = 0;

   D = new tensor**[SIZE];

   for(i = 0; i < SIZE; i++)
     D[i] = new tensor*[SIZE];

   for(i = 0; i < SIZE; i++)
     for(j = 0; j < SIZE; j++)
       D[i][j] = new tensor[SIZE];

   for(i = 0; i < SIZE; i++)
     for(j = 0; j < SIZE; j++)
       for(k = 0; k < SIZE; k++)
         D[i][j][k].Resize(2,2,2);



}
