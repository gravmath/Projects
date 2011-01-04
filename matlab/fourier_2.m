function fourier_2(N,lim,step)

x = lim(1):step:lim(2);


f = zeros(1,max(size(x)));

for j = 1:1:N
   expon    = -1i*pi*j;
   a        =  1i*pi*j;
   coeff_1  = 2/(a^2)*( 1 - (a+1)*exp(expon) );
   coeff_2  = (exp(2*expon)-exp(expon))/a;
   coeff_np = (coeff_1 + coeff_2)/2;
   f = f + coeff_np*exp(1i*pi*j*x);
end
for j = -N:1:-1
   expon    = -1i*pi*j;
   a        =  1i*pi*j;
   coeff_1  = 2/(a^2)*( 1 - (a+1)*exp(expon) );
   coeff_2  = (exp(2*expon)-exp(expon))/a;
   coeff_np = (coeff_1 + coeff_2)/2;
   f = f + coeff_np*exp(1i*pi*j*x);
end
   
plot(x,f,'bo-')