function dS = MSP(t,S,parms)

%unpack physical parameters
mu = parms.mu;
J2 = parms.J2;
R  = parms.R;

%allocate space 
dS        = zeros(6,1);

%equate the time derivative of the position with the velocity
dS(1:3,1) = S(4:6,1);     

%unpack the postion to temporary variables for ease of reading/maintenance
x         = S(1);
y         = S(2);
z         = S(3);
r         = sqrt( x*x + y*y + z*z );

dS(4,1)   = -mu*x/r^3 * (1 - J2*3/2*(R/r)^2*(5*z^2/r^2 - 1));
dS(5,1)   = y/x*dS(4);
dS(6,1)   = -mu*z/r^3 * (1 + J2*3/2*(R/r)^2*(3 - 5*z^2/r^2));

