close all;
options = odeset('RelTol',1e-11);

parms.k  = 2;
parms.L0 = 6.8;
parms.m  = 1/2;
parms.h  = 7;

linear_osc_freq  = sqrt(2*parms.k/parms.m);
first_order_freq = linear_osc_freq*sqrt(1-parms.L0/parms.h);
period           = 2*pi/linear_osc_freq;

time_span = [0:0.1:10*period];
x0        = 2;
v0        = 0;
ICs       = [x0,v0];

[T,S]                  = ode45('slidder_ode',time_span,ICs,options,parms); 
lin_osc_x              = ICs(1)*cos(linear_osc_freq*T)  + ICs(2)*sin(linear_osc_freq*T);  
first_ord_correction_x = ICs(1)*cos(first_order_freq*T) + ICs(2)*sin(first_order_freq*T);

plot(T,S(:,1),'ko-',T,lin_osc_x,'r-',T,first_ord_correction_x,'b-')
legend('Numerical','Linear Osc.','1st-order corr.');