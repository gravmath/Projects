function fourier_4(N,lim,step)

x = lim(1):step:lim(2);


f = zeros(1,max(size(x)))+0.5;

for j = 1:1:N
   m        = j;
   coeff_p  = ((2i*m*pi+1)*exp(-2i*m*pi) - 1)/(4*m^2*pi^2);
   m        = -j;
   coeff_m  = ((2i*m*pi+1)*exp(-2i*m*pi) - 1)/(4*m^2*pi^2);
   expon_p  =  1i*2*pi*j;
   expon_m  = -1i*2*pi*j;
      
   f = f + coeff_p*exp(expon_p*x) + coeff_m*exp(expon_m*x);
end
   
plot(x,f,'bo-')