function ds = slidder_ode(t,s,dummy,parms)

k  = parms.k;
L0 = parms.L0;
m  = parms.m;
h  = parms.h;

L  = sqrt( s(1)^2 + h^2 );

ds(1) = s(2);
ds(2) = -2*k/m*s(1)*(1-L0/L);
ds    = ds';


