#include <math.h>
#define STD_SCHW 1
#define ISO_SCHW 2

void     main(void);
void       std_schw(double b[3], double *a,     double  g[3][3]);
void     std_schw_d(double b[3], double da[3], double dg[3][3][3]);
void       iso_schw(double b[3], double *a,     double  g[3][3]);
void     iso_schw_d(double b[3], double da[3], double dg[3][3][3]);
double   W(double a[3], double b[3], double h);
void     dW(double a[3], double b[3], double h, double ret[3]);


void     main(void)
{
   int metric_choice;
   metric_choice = STD_SCHW;
   double b[3], min[3], delta, h;
   double z[3], u[3],   dr[3], Du[3], du[3], dstate[6], Dstate[6];
   int num;
   z[0]   = 20.0;
   z[1]   = 0.0;
   z[2]   = 0.0;
   u[0]   = 0.000001533;
   u[1]   = 0.242551449;
   u[2]   = 0.0;
   h      = 1.0/32.0;
   min[0] = -21.199999999999999000;
   min[1] = -21.199999999999999000;
   min[2] = -1.231535269709543400;
   delta  = 0.175933609958506210;

   min[0] = z[0] - 1.0*h;
   min[1] = z[1] - 1.0*h;
   min[2] = z[2] - 1.0*h;
   num    = 200;
   delta  = 2.0*h/( (double) num );
   double w3,      dw3[3],     sum3,           dsum3[3];
   double a,       a_sm,       da[3],          Da_sm[3],       da_sm[3];
   double g[3][3], g_sm[3][3], dg[3][3][3],    Dg_sm[3][3][3], dg_sm[3][3][3];
   double detg,    r3g;
   double detG,    Ginv[3][3], DGinv[3][3][3], dGinv[3][3][3];
   double Lambda,  inv_L,      H,              dummy,          dummy_v[3];


   int    i, j, k;            //spatial indices
   int    m, p, q, r, s, t;   //tensor indices

   sum3     = 0.0;
   a_sm     = 0.0;
   for(m=0;m<3;m++)
     {
       dsum3[m] = 0.0;
       Da_sm[m] = 0.0;
       da_sm[m] = 0.0;
     }
    for(m=0;m<3;m++) for(p=0;p<3;p++) g_sm[m][p] = 0.0;
    for(m=0;m<3;m++) for(p=0;p<3;p++) for(q=0;q<3;q++)
    {
       Dg_sm[m][p][q] = 0.0;
       dg_sm[m][p][q] = 0.0;
    }

   for(i=0; i < num; i++)
     for(j=0; j < num; j++)
       for(k=0; k < num; k++)
         {
            b[0] = min[0] + delta * (double)i;
            b[1] = min[1] + delta * (double)j;
            b[2] = min[2] + delta * (double)k;
            w3 = W(z,b,1.0);
            if( w3 > 0.0 )
            {
               dW(z,b,1.0,dw3);
               switch (metric_choice)
               {
                 case STD_SCHW:
                   std_schw(b,&a,g);
                   std_schw_d(b,da,dg);
                   break;
                 case ISO_SCHW:
                   iso_schw(b,&a,g);
                   iso_schw_d(b,da,dg);
                   break;
               }
               detg =   g[0][0]*(g[1][1]*g[2][2] - g[1][2]*g[2][1])
                      + g[0][1]*(g[2][0]*g[1][2] - g[1][0]*g[2][2])
                      + g[0][2]*(g[1][0]*g[2][1] - g[1][1]*g[2][0]);
               r3g  = sqrt(detg);
               r3g  = 1.0;

               sum3 = sum3 + w3*r3g;
               a_sm = a_sm + a*w3*r3g;
               for(m=0;m<3;m++) dsum3[m] += dw3[m];
               for(m=0;m<3;m++) Da_sm[m] += da[m]*w3*r3g;
               for(m=0;m<3;m++) da_sm[m] += a*dw3[m];
               for(m=0;m<3;m++) for(p=0;p<3;p++)  g_sm[m][p] += g[m][p]*w3*r3g;
               for(m=0;m<3;m++) for(p=0;p<3;p++) for(q=0;q<3;q++)
                   Dg_sm[m][p][q] += dg[m][p][q]* w3 * r3g;
               for(m=0;m<3;m++) for(p=0;p<3;p++) for(q=0;q<3;q++)
                   dg_sm[m][p][q] += g[m][p]*dw3[q];
            }
         }
a_sm = a_sm / sum3;
for(m=0;m<3;m++) Da_sm[m] = Da_sm[m]/sum3;
for(m=0;m<3;m++) da_sm[m] = da_sm[m]/sum3;
for(m=0;m<3;m++) da_sm[m] = da_sm[m] - a_sm * dsum3[m] / sum3;
for(m=0;m<3;m++) for(p=0;p<3;p++) g_sm[m][p] = g_sm[m][p]/sum3;
for(m=0;m<3;m++) for(p=0;p<3;p++) for(q=0;q<3;q++)
  Dg_sm[m][p][q] = Dg_sm[m][p][q]/sum3;
for(m=0;m<3;m++) for(p=0;p<3;p++) for(q=0;q<3;q++)
  dg_sm[m][p][q] = dg_sm[m][p][q]/sum3 - g_sm[m][p]*dsum3[q] / sum3;

 //Form the inverse of the smoothed metric
  detG =   g_sm[0][0]*(g_sm[1][1]*g_sm[2][2] - g_sm[1][2]*g_sm[2][1])
         + g_sm[0][1]*(g_sm[2][0]*g_sm[1][2] - g_sm[1][0]*g_sm[2][2])
         + g_sm[0][2]*(g_sm[1][0]*g_sm[2][1] - g_sm[1][1]*g_sm[2][0]);

  Ginv[0][0] = ( g_sm[1][1]*g_sm[2][2] - g_sm[1][2]*g_sm[2][1] )/detG;
  Ginv[0][1] = ( g_sm[0][2]*g_sm[2][1] - g_sm[0][1]*g_sm[2][2] )/detG;
  Ginv[0][2] = ( g_sm[0][1]*g_sm[1][2] - g_sm[0][2]*g_sm[1][1] )/detG;
  Ginv[1][0] = ( g_sm[1][2]*g_sm[2][0] - g_sm[1][0]*g_sm[2][2] )/detG;
  Ginv[1][1] = ( g_sm[0][0]*g_sm[2][2] - g_sm[0][2]*g_sm[2][0] )/detG;
  Ginv[1][2] = ( g_sm[0][2]*g_sm[1][0] - g_sm[0][0]*g_sm[1][2] )/detG;
  Ginv[2][0] = ( g_sm[1][0]*g_sm[2][1] - g_sm[1][1]*g_sm[2][0] )/detG;
  Ginv[2][1] = ( g_sm[0][1]*g_sm[2][0] - g_sm[0][0]*g_sm[2][1] )/detG;
  Ginv[2][2] = ( g_sm[0][0]*g_sm[1][1] - g_sm[0][1]*g_sm[1][0] )/detG;


  //now calculate the tensor needed for the equations of motion
  for(m=0;m<3;m++)  for(p=0;p<3;p++)  for(q=0;q<3;q++)
  {
     DGinv[m][p][q] = 0.0;
     dGinv[m][p][q] = 0.0;
     for(r=0;r<3;r++) for(s=0;s<3;s++)
     {
        DGinv[m][p][q] += -1.0 * Ginv[m][r] * Dg_sm[r][s][q] * Ginv[s][p];
        dGinv[m][p][q] += -1.0 * Ginv[m][r] * dg_sm[r][s][q] * Ginv[s][p];
     }
  }

  //now calculate the RHS
  dummy = 0.0;
  for(m=0;m<3;m++) for(p=0;p<3;p++) dummy += u[m] * Ginv[m][p] * u[p];
  Lambda = sqrt(1.0 + dummy);
  inv_L = 1.0/Lambda;

  //dr
  for(m=0;m<3;m++)
  {
    dr[m] = 0.0;
    for(p=0;p<3;p++) dr[m] += a_sm * inv_L * Ginv[m][p] * u[p];
  }

  //du and Du
  for(m=0;m<3;m++)
  {
     dummy_v[m] = 0.0;
     for(p=0;p<3;p++) for(q=0;q<3;q++)
        dummy_v[m] += u[p] * dGinv[p][q][m] * u[q];
     du[m] = -Lambda * da_sm[m] - a_sm * 0.5 * inv_L *dummy_v[m];

     dummy_v[m] = 0.0;
     for(p=0;p<3;p++) for(q=0;q<3;q++)
        dummy_v[m] += u[p] * DGinv[p][q][m] * u[q];
     Du[m] = -Lambda * da_sm[m] - a_sm * 0.5 * inv_L *dummy_v[m];
  }

  //Calculate the conserved quantities
  H =  a_sm*Lambda;

  //pack the state
  for( i = 0; i < 3; i++)
  {
    dstate[i]   = dr[i];
    dstate[i+3] = du[i];
    Dstate[i]   = dr[i];
    Dstate[i+3]   = Du[i];
  }


}


//**********************************************************************
//std_schw
//**********************************************************************
void std_schw(double b[3], double *a, double g[3][3])
{
   double M, r, M2_r, Q;

   M = 1.0;
   r = sqrt( b[0]*b[0] + b[1]*b[1] + b[2]*b[2] );
   M2_r = 2.0*M/r;
   Q = M2_r/( r * r * ( 1.0 - M2_r ) );

   *a = sqrt(1.0 - M2_r);

   g[0][0] = 1.0 + b[0]*b[0]*Q;
   g[0][1] =       b[0]*b[1]*Q;
   g[0][2] =       b[0]*b[2]*Q;
   g[1][0] =       b[1]*b[0]*Q;
   g[1][1] = 1.0 + b[1]*b[1]*Q;
   g[1][2] =       b[1]*b[2]*Q;
   g[2][0] =       b[2]*b[0]*Q;
   g[2][1] =       b[2]*b[1]*Q;
   g[2][2] = 1.0 + b[2]*b[2]*Q;
}

//**********************************************************************
//std_schw_d
//**********************************************************************
void std_schw_d(double b[3], double da[3], double dg[3][3][3])
{
   double M2, r, x, y, z, inv_r, inv_r3, a, inv_a, Q, P;
   int i,j,k;

   M2     = 2.0;
   r      = sqrt( b[0]*b[0] + b[1]*b[1] + b[2]*b[2] );
   inv_r  = 1.0/r;
   inv_r3 = inv_r * inv_r * inv_r;
   Q      = M2*inv_r3/(1.0 - M2*inv_r);
   P      = (3.0*inv_r*inv_r + Q);
   x      = b[0];
   y      = b[1];
   z      = b[2];

   dg[0][0][0] = 2.0*x - P*x*x*x;
   dg[0][0][1] =       - P*x*x*y;
   dg[0][0][2] =       - P*x*x*z;
   dg[0][1][0] =     y - P*x*x*y;
   dg[0][1][1] =     x - P*x*y*y;
   dg[0][1][2] =       - P*x*y*z;
   dg[0][2][0] =     z - P*x*x*z;
   dg[0][2][1] =       - P*x*y*z;
   dg[0][2][2] =     x - P*x*z*z;

   dg[1][0][0] =     y - P*x*x*y;
   dg[1][0][1] =     x - P*x*y*y;
   dg[1][0][2] =       - P*x*y*z;
   dg[1][1][0] =       - P*x*y*y;
   dg[1][1][1] = 2.0*y - P*y*y*y;
   dg[1][1][2] =       - P*y*y*z;
   dg[1][2][0] =       - P*x*y*z;
   dg[1][2][1] =     z - P*y*y*z;
   dg[1][2][2] =     y - P*y*z*z;

   dg[2][0][0] =     z - P*x*x*z;
   dg[2][0][1] =       - P*x*y*z;
   dg[2][0][2] =     x - P*x*z*z;
   dg[2][1][0] =       - P*x*x*z;
   dg[2][1][1] =     z - P*y*y*z;
   dg[2][1][2] =     y - P*y*z*z;
   dg[2][2][0] =       - P*x*z*z;
   dg[2][2][1] =       - P*y*z*z;
   dg[2][2][2] = 2.0*z - P*z*z*z;

   for(i=0;i<3;i++) for(j=0;j<3;j++) for(k=0;k<3;k++) dg[i][j][k] *= Q;
   a     = sqrt(1.0 - M2*inv_r);
   inv_a = 1.0/a;
   da[0] = 0.5*M2*inv_r3*inv_a*x;
   da[1] = 0.5*M2*inv_r3*inv_a*y;
   da[2] = 0.5*M2*inv_r3*inv_a*z;
}

//**********************************************************************
//iso_schw
//**********************************************************************
void iso_schw(double b[3], double *a, double g[3][3])
{
   double M, r, M_2r, op_M_2r, om_M_2r, Q, R;

   M = 1.0;
   r = sqrt( b[0] * b[0] + b[1] * b[1] + b[2] * b[2] );
   M_2r    = 0.5*M/r;
   op_M_2r = 1.0 + M_2r;
   om_M_2r = 1.0 - M_2r;
   Q       = om_M_2r / op_M_2r;
   R       = op_M_2r * op_M_2r * op_M_2r * op_M_2r;

   *a = Q;

   g[0][0] = R;
   g[0][1] = 0.0;
   g[0][2] = 0.0;
   g[1][0] = 0.0;
   g[1][1] = R;
   g[1][2] = 0.0;
   g[2][0] = 0.0;
   g[2][1] = 0.0;
   g[2][2] = R;

}

//**********************************************************************
//iso_schw_d
//**********************************************************************
void iso_schw_d(double b[3], double da[3], double dg[3][3][3])
{
   double M, M_2r, r, x, y, z, inv_r, inv_r3, op_M_2r, inv_op_2, M_inv_r3, op_3;

   M = 1.0;
   r = sqrt( b[0] * b[0] + b[1] * b[1] + b[2] * b[2] );

   M_2r     = 0.5*M/r;
   inv_r    = 1.0/r;
   inv_r3   = inv_r * inv_r * inv_r;
   op_M_2r  = 1.0 + M_2r;
   op_3     = op_M_2r * op_M_2r * op_M_2r;
   inv_op_2 = 1.0/op_M_2r/op_M_2r;
   M_inv_r3 = M*inv_r3;
   x = b[0];
   y = b[1];
   z = b[2];

   da[0] = x*M_inv_r3*inv_op_2;
   da[1] = y*M_inv_r3*inv_op_2;
   da[2] = z*M_inv_r3*inv_op_2;

   dg[0][0][0] = -2.0*op_3*M_inv_r3*x;
   dg[0][0][1] = -2.0*op_3*M_inv_r3*y;
   dg[0][0][2] = -2.0*op_3*M_inv_r3*z;
   dg[0][1][0] =  0.0;
   dg[0][1][1] =  0.0;
   dg[0][1][2] =  0.0;
   dg[0][2][0] =  0.0;
   dg[0][2][1] =  0.0;
   dg[0][2][2] =  0.0;

   dg[1][0][0] =  0.0;
   dg[1][0][1] =  0.0;
   dg[1][0][2] =  0.0;
   dg[1][1][0] = -2.0*op_3*M_inv_r3*x;
   dg[1][1][1] = -2.0*op_3*M_inv_r3*y;
   dg[1][1][2] = -2.0*op_3*M_inv_r3*z;
   dg[1][2][0] =  0.0;
   dg[1][2][1] =  0.0;
   dg[1][2][2] =  0.0;

   dg[2][0][0] =  0.0;
   dg[2][0][1] =  0.0;
   dg[2][0][2] =  0.0;
   dg[2][1][0] =  0.0;
   dg[2][1][1] =  0.0;
   dg[2][1][2] =  0.0;
   dg[2][2][0] = -2.0*op_3*M_inv_r3*x;
   dg[2][2][1] = -2.0*op_3*M_inv_r3*y;
   dg[2][2][2] = -2.0*op_3*M_inv_r3*z;
}

//**********************************************************************
//W
//**********************************************************************
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

//**********************************************************************
//dW
//**********************************************************************
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

