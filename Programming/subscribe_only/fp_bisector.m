%fp_bisector - designed to run fp_runner and to target on R = constant and sigmaR < 1.5e-4
function bisect_flag = fp_bisector(std_radius,metric,kernel,h,filename)

P1 = path;
path(P1,'\school\analysis\matlab\');
M             = 1;
u             = 0;
kernel_choice = kernel;
num_iters     = 30;
num_pts       = 2000;
num_orbits    = 1;

%create first guess for v
if strcmp(metric,'STD_SCHW')
   radius = std_radius;
   v      = 1/sqrt(radius - 3);
elseif strcmp(metric,'ISO_SCHW')
   radius = ( (std_radius - M) + sqrt( std_radius * (std_radius - 2*M) ) )/2;
   v      = sqrt( (2*radius + 1)^4/radius/(4*radius^2 -8*radius+1))/2/radius;
end

delta  = 1e-2;

%Now run the first 'A' trial
init_A   = 0;
counterA = 1;
deltaA   = delta;
while( init_A == 0 & counterA < 30 )
   vA = v - deltaA;
   [RA,SA,minA,maxA] = fp_runner(std_radius,metric,kernel,h,num_pts,'r3g',num_orbits,filename,u,vA,'KD');
   dist_min = radius - minA;
   dist_max = radius - maxA;
   dist_RA  = radius - RA;
   if( dist_RA  > 1e-4 )
      init_A = 1; %velocity bracketed
   elseif( minA < radius & maxA < radius ) %trending down
      init_A = 1;
   else
      deltaA = 2.0*deltaA;
   end
   counterA = counterA + 1;
end

%Now run the first 'B' trial
init_B   = 0;
counterB = 1;
deltaB   = delta;
while( init_B == 0 & counterB < 30 )
   vB = vA + deltaB/2;
   [RB,SB,minB,maxB] = fp_runner(std_radius,metric,kernel,h,num_pts,'r3g',num_orbits,filename,u,vB,'KD');
   dist_min = radius - minB;
   dist_max = radius - maxB;
   dist_RB  = radius - RB;
   if( dist_RB  < -1e-4 )
      init_B = 1; %velocity bracketed
   elseif( minB > radius & maxB > radius ) %trending up
      init_B = 1;
   else
      deltaB = 1.1*deltaB;
   end
   counterB = counterB + 1;
end

%Finally bisect
bisect_flag = -1;
continue    = 1;
ST          = 1;
if( init_A == 1 & init_B == 1 )
   bisect_flag = 0;
   i = 0;
   while( continue == 1 )
      vT = (vB + vA)/2;
      [RT,ST,minT,maxT] = fp_runner(std_radius,metric,kernel,h,num_pts,'r3g',num_orbits,filename,u,vT,'KD');
      i;
      vA;
      vB;
      ST;
      dist_min = radius - minT;
      dist_max = radius - maxT;
      dist_RT  = radius - RT;
      i = i + 1;
      if( ST  < 1e-5 )
         continue = 0; %all done
      elseif( dist_min <= 0 & dist_max <= 0 ) %T is above and velocity bracketed by A and T
         vB = vT;
      elseif( dist_min >= 0 & dist_max >= 0 ) %T is below and velocity bracketed by T and B
         vA = vT;
      elseif( dist_min >= 0 & dist_max <= 0 ) %T straddles the radius 
         if(  RT > radius )
            vB = vT;
         else
            vA = vT;
         end
      end
      if( i > num_iters )
         continue = 0;
      end
   end
end
if( i - 1 < num_iters ) 
   bisect_flag = 1;
end
bisect_flag
%Now produce the 10 orbit sample if the bisection was successful
num_orbits = 10;
[R,S,min,max] = fp_runner(std_radius,metric,kernel,h,num_pts,'r3g',num_orbits,filename,u,vT,'KD');

   
   
         
      

