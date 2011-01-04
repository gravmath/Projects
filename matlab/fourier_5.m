function fourier_5(N,lim,step)

x = lim(1):step:lim(2);

v = 1/6;
A = -2i*pi*v;

f = zeros(1,max(size(x))) + 36*v;
for j = 1:1:N
   m        = j;
   coeff_p   = exp(-3*m*A)*( (6*m*A - 2)*exp(6*m*A) + 6*m*A + 2 )/(m^3*A^3)*v;
   m        = -j;
   coeff_m   = exp(-3*m*A)*( (6*m*A - 2)*exp(6*m*A) + 6*m*A + 2 )/(m^3*A^3)*v;
   f = f + coeff_p*exp(j*A*x) + coeff_m*exp(-j*A*x);
end
   
plot(x,f,'bo-')