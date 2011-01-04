function fourier_1(N,lim,step)

x = lim(1):step:lim(2);

a = 0.5;
b = 1.0;

%N = 1;

f = ones(1,max(size(x)))*2*a/(a+b);

for j = 1:1:N
   coeff_n = sin(2*pi*j*a/(a+b))/pi/j;
   f = f + coeff_n*exp(2i*pi*j*x/(a+b)) + coeff_n*exp(-2i*pi*j*x/(a+b));
end
plot(x,f)