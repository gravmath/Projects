//Tensor and array test - version 2
//
//Purpose:	Test the concept of making a tensor field
//
//Date:		4/7/03


#include "C:\My Documents\Programming\tensor++\tensor.hpp"

void main(void);

void main(void)
{
   tensor skippy(2,3,3);
   tensor *dude;
   char *dump;
   FILE *file;
   dump = NULL;

   //Just stretching my legs here
   for(int i = 0; i < 3; i++)
     for(int j = 0; j < 3; j++)
       skippy.Set(i*i + j,i,j);

   skippy.print("skippy",&dump);
   file = fopen("test2.txt","w");
   fprintf(file,"The formula for skippy is skippy(i,j) = i*i +j\n");
   fprintf(file,"%s",dump);
   fprintf(file,"\n",dump);

   //Define a 1-d field of 2x2 matrices
   dude = new tensor[5];
   dude[0].Resize(2,2,2);
   dude[1].Resize(2,2,2);
   dude[2].Resize(2,2,2);
   dude[3].Resize(2,2,2);
   dude[4].Resize(2,2,2);

   for(int i = 0; i < 5; i++)
     for(int j = 0; j < 2; j++)
       for(int k = 0; k < 2; k++)
         dude[i].Set(i+j+k,j,k);

   fprintf(file,
   "The formula for dude (a 1-D field of 2x2 matrices) is dude[i](j,k) = i+j+k\n");
            
   for(int i = 0; i < 5; i++)
	{
     dude[i].print("dude[i]",&dump);
     fprintf(file,"i = %d\n",i);
     fprintf(file,dump);
   }
   fprintf(file,"\n");

   tensor test;
   fprintf(file,"The formula for test is test = dude[i] + dude[i+1] i = 0..3\n");
   for(int i = 0; i < 4; i++)
   {
     test <= dude[i] + dude[i+1];
     fprintf(file,"i = %d\n",i);     
     test.print("test",&dump);
     fprintf(file,"%s",dump);
   }
}
