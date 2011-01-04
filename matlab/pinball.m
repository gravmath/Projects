clear all
close all

N      = 20;
num_MC = 5000;
p_r    = 29/59;
p_d    = 1/59;
p_l    = 29/59;

%p_r    = 1/3;
%p_d    = 1/3;
%p_l    = 1/3;

x = zeros(2*N+1,1);
for j = 1:1:num_MC
   current_loc = 0;
   for i = 1:1:N
      prob = rand(1,1);
      if( prob <= p_r )
         current_loc = current_loc + 1;
      elseif( prob > p_r & prob <= (p_r + p_d) )
         current_loc = current_loc;
      else
         current_loc = current_loc - 1;
      end
   end
   x(current_loc + N + 1) = x(current_loc + N + 1) + 1;
end
x = x/num_MC;
m = -N:1:N;

plot(m,x)
   
         
      
   