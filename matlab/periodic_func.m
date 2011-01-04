%periodic_func

function periodic_func(min_q,max_q)

q          = min_q:0.01:max_q;
saw_period = 2;
saw_offset = -1;

for i = 1:1:max(size(q))
   s(i) = saw_tooth(q(i),saw_period,saw_offset);
end

plot(q,s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ret_val = saw_tooth(x,period,offset)

x = mod(x-offset,period);

if( x >= 0 & x <= period/2 )
   ret_val = x;
elseif( x >= period/2 & x <= period )
   ret_val = period - x;
end
