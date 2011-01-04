function fourier_3(N,lim,step)

x = lim(1):step:lim(2);


f = zeros(1,max(size(x)))+2;

for j = 1:1:N
   coeff    = -9/4/pi^2/j^2;
   expon_p  =  1i*2*pi*j/3;
   expon_m  = -1i*2*pi*j/3;
   
   coeff_p  = coeff*(1 - exp(expon_m) - exp(2*expon_m) + exp(3*expon_m) );
   coeff_m  = coeff*(1 - exp(expon_p) - exp(2*expon_p) + exp(3*expon_p) );
      
   f = f + coeff_p*exp(expon_p*x) + coeff_m*exp(expon_m*x);
end
   
plot(x,f,'bo-')