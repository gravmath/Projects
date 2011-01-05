#include <math.h>

void     main(void);
double   F(double b[3]);
void     dF(double b[3],double ret[3]);
double   W(double a[3], double b[3], double h);
void     dW(double a[3], double b[3], double h, double ret[3]);


void     main(void)
{
   double a[3], min[3], delta, h;
   int num;
   a[0]   = 20.0;
   a[1]   = 0.0;
   a[2]   = 0.0;
   h      = 1.0;
   //min[0] = -21.199999999999999000;
   //min[1] = -21.199999999999999000;
   //min[2] = -1.231535269709543400;
   //delta  = 0.175933609958506210;

   min[0] = a[0] - 1.0*h;
   min[1] = a[1] - 1.0*h;
   min[2] = a[2] - 1.0*h;
   delta  = 0.01;
   num    = 2.0*h/delta;
   double sum3, dsum3[3], f3, Df3[3], df3[3], b[3], w3, f, df[3], dw3[3];
   int    i,    j,        k,  m;

   sum3     = 0.0;
   for(m=0;m<3;m++)
     {
       dsum3[m] = 0.0;
       Df3[m]   = 0.0;
       df3[m]   = 0.0;
     }
   f3       = 0.0;

   for(i=0; i < num; i++)
     for(j=0; j < num; j++)
       for(k=0; k < num; k++)
         {
            b[0] = min[0] + delta * (double)i;
            b[1] = min[1] + delta * (double)j;
            b[2] = min[2] + delta * (double)k;
            w3 = W(a,b,1.0);
            if( w3 > 0.0 )
            {
               dW(a,b,1.0,dw3);
               f    =  F(b);
               dF(b,df);
               sum3 = sum3 + w3;
               for(m=0;m<3;m++) dsum3[m] = dsum3[m] + dw3[m];
               f3   = f3 + w3*f;
               for(m=0;m<3;m++) Df3[m]   = Df3[m] + df[m]*w3;
               for(m=0;m<3;m++) df3[m]   = df3[m] + f*dw3[m];
            }
         }
f3 = f3/sum3;
for(m=0;m<3;m++) Df3[m] = Df3[m]/sum3;
for(m=0;m<3;m++) df3[m] = df3[m]/sum3;
for(m=0;m<3;m++) df3[m] = df3[m] - f3 * dsum3[m] / sum3;
}

double   F(double b[3])
{
   double r;
   r = sqrt( b[0]*b[0] + b[1]*b[1] + b[2]*b[2] );
   return sqrt(1.0 - 2.0 /r);

}

void dF(double b[3], double ret[3])
{
   double r, a;
   r = sqrt( b[0]*b[0] + b[1]*b[1] + b[2]*b[2] );
   a = sqrt(1.0 - 2.0 /r);
   ret[0] = b[0] / a / r / r / r;
   ret[1] = b[1] / a / r / r / r;
   ret[2] = b[2] / a / r / r / r;
}

double   W(double a[3], double b[3], double h)
{
   double v[3], r, pi;

   v[0] = a[0] - b[0];
   v[1] = a[1] - b[1];
   v[2] = a[2] - b[2];

   r = sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );

   pi = 3.14159265358979;

   if( r > h )
     return 0.0;
   else
     return 315.0/64.0/pi/h/h/h *
              (1.0 - r*r/h/h) * (1.0 - r*r/h/h) * (1.0 - r*r/h/h);
}

void dW(double a[3], double b[3], double h, double ret[3])
{
   double v[3], r, pi;

   v[0] = a[0] - b[0];
   v[1] = a[1] - b[1];
   v[2] = a[2] - b[2];

   r = sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );

   pi = 3.14159265358979;

   if( r > h )
   {
      ret[0] = 0.0;
      ret[1] = 0.0;
      ret[2] = 0.0;
   }
   else
   {
      ret[0] = -945.0/32.0/pi/h/h/h/h/h *
                (1.0 - r*r/h/h) * (1.0 - r*r/h/h) * v[0];
      ret[1] = -945.0/32.0/pi/h/h/h/h/h *
                (1.0 - r*r/h/h) * (1.0 - r*r/h/h) * v[1];
      ret[2] = -945.0/32.0/pi/h/h/h/h/h *
                (1.0 - r*r/h/h) * (1.0 - r*r/h/h) * v[2];
   }
}

