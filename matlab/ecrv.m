%initialize
clear all
close all

%define the parameters of the sim
q          = 1;
p          = 1;
tau        = 100;
dt         = 1;
num_trials = 5000;
num_steps  = 300;
offset     = 1;
time       = 0:dt:(num_steps-1)*dt;

%define exponential coefficients
w1 = exp(-dt/tau);
w2 = sqrt( q*tau/2 * (1 - exp(-2*dt/tau)));

%allocate space for the positions
x = zeros(num_trials,num_steps);

%define the initial positions
x(:,1) = -3*p + 6*p*rand(num_trials,1) + offset;

subplot(2,1,1)
for i = 1:num_trials
    for j = 2:num_steps
        x(i,j) = w1*x(i,j-1) + w2*randn(1);
    end
    hold on
    plot(time,x(i,:),'Color',[rand(1),rand(1),rand(1)]);
end
hold off

mean_x0 = mean(x(:,1));
p0      = std(x(:,1))^2;

subplot(2,2,3);
plot(time,mean(x),time,mean_x0*exp(-time/tau),'r-');

subplot(2,2,4);
plot(time,std(x),time,sqrt(exp(-2*time/tau)*p0+q*tau/2*(1-exp(-2*time/tau))),'r-');
