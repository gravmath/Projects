close all;
clear all;

options = odeset('RelTol',1e-11);

parms.mu  = 398600.4414;
parms.R   = 6378.14;
parms.J2  = 1082.64e-6;
parms.J2  = 0;

min       = 60;
time_span = [0:min:100*min];
ICs       = [7000;0;0;0;7.54/sqrt(2);7.54/sqrt(2)];

[T,S]     = ode45('MSP',time_span,ICs,options,parms); 

