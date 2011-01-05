function [t_star,energy] = find_period(A)

a = A(1,3);
for i = 1:1:max(size(A))
   b = A(i,3);
   sab = sign( a * b);
   if( sab < 0 )
      if( b > 0 )
         counter = i;
         break;
      end
   end
   a = b;
end
time   = A(counter-10:1:counter+10,1);
y      = A(counter-10:1:counter+10,3);
t_star = spline(y,time,0);
energy = mean(A(2:counter,8));