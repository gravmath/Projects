%fp_runner
function [Rbar,Rsigma,Rmin,Rmax] = fp_runner(std_radius,metric,kernel,h,num_pts,slicing,num_period,filename,u,v,derivs_choice)

P1 = path;
path(P1,'c:\school\analysis\matlab\');
M             = 1;
kernel_choice = kernel;

%first create the first guess
Period = 2*pi*sqrt(std_radius^3/M);
Period = ceil(Period);

if strcmp(metric,'STD_SCHW')
   radius = std_radius;
   if( v == -1 )
      v      = 1/sqrt(radius - 3);
   end
elseif strcmp(metric,'ISO_SCHW')
   radius = ( (std_radius - M) + sqrt( std_radius * (std_radius - 2*M) ) )/2;
   if( v == -1 )
      v      = sqrt( (2*radius + 1)^4/radius/(4*radius^2 -8*radius+1))/2/radius;
   end
end
[min_x,min_y,min_z,nx,ny,nz,delta] = set_sim_parms(radius,h,num_pts);

time(1)  = 0.1;
time(3)  = 50;
N        = ceil( ceil(Period)/time(1)/time(3) * num_period ) + 1;
time(2)  = (N - 1) * time(3);
vel      = [u;v];
grid(1)  = min_x;
grid(2)  = min_y;
grid(3)  = min_z;
grid(4)  = nx;
grid(5)  = ny;
grid(6)  = nz;
grid(7)  = delta;
grid(8)  = delta;
grid(9)  = delta;
grid(10) = num_pts;

write_FP_file('C:\School\Programming\subscribe_only\subscribe_geodesics_parms.txt',M,radius,metric,kernel_choice,h,time,vel,grid,derivs_choice);
%run nominal
if strcmp(slicing,'bare')
   !subscribe_geodesics_bare
elseif strcmp(slicing,'r3g')
   !subscribe_geodesics_r3g
elseif strcmp(slicing,'spherical')
   !subscribe_geodesics_spherical
elseif strcmp(slicing,'scalar')
   !subscribe_geodesics_scalar
end   
clear ephem;
clear R;
ephem      = load('C:\school\programming\subscribe_only\ephem.txt');
time       = ephem(:,1);
R          = sqrt( ephem(:,2).*ephem(:,2) + ephem(:,3).*ephem(:,3) + ephem(:,4).*ephem(:,4) );
Rbar       = mean(R);
Rsigma     = std(R);
%Rhalf      = R(floor(max(size(R)))/2); %Radius value at the half orbit point
Rmin       = min(R);
Rmax       = max(R);
plot(time,R)
string = ['Subscribe Geo Run: R =',num2str(radius),' Rbar and Rsigma = ',num2str(Rbar),' and ',num2str(Rsigma),...
          ' h = ',num2str(h),' K = ',num2str(kernel_choice),' metric = ',metric,' Vel = ',num2str(u),' and ',num2str(v)];
title(string);
dos(['copy ephem.txt e:\kd_runs\',filename]);
fig_name = filename(1:max(size(filename))-4);
fig_name = ['e:\kd_runs\',fig_name,'.jpg'];
print('-djpeg',fig_name);