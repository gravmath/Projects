%fp_scanner
function fp_scanner(std_radius,metric,kernel,h,num_pts,slicing,max_iters,filename)

P1 = path;
path(P1,'\school\analysis\matlab\');
fid = fopen('targ_hist','a');
M             = 1;
u             = 0;
kernel_choice = kernel;
du            = 1e-4;
dv            = 1e-4;
num_iters     = 10;
Rbar          = 1;
Rsigma        = 1;
tolRbar       = 1e-6;
tolRsigma     = 1e-6;
best_u        = 0;
best_v        = 0;
best_Rbar     = 100;
best_Rsigma   = 100;

%first create the first guess
if strcmp(metric,'STD_SCHW')
   radius = std_radius;
   v      = 1/sqrt(radius - 3);
elseif strcmp(metric,'ISO_SCHW')
   radius = ( (std_radius - M) + sqrt( std_radius * (std_radius - 2*M) ) )/2;
   v      = sqrt( (2*radius + 1)^4/radius/(4*radius^2 -8*radius+1))/2/radius;
end
[min_x,min_y,min_z,nx,ny,nz,delta] = set_sim_parms(radius,h,num_pts);

time(1)  = 0.1;
time(2)  = 10000;
time(3)  = 50;
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

counter = 1;
diff_old = [1;1];
while ( counter  < max_iters )
   counter
   vel(1) = u % - 0.0001 + (counter - 1)*0.00004;
   vel(2) = v - 0.0001 + (counter - 1)*0.00004
   write_FP_file('\School\Programming\subscribe_only\subscribe_geodesics_parms.txt',M,radius,metric,kernel_choice,h,time,vel,grid);
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
   ephem      = load('\school\programming\subscribe_only\ephem.txt');
   R          = sqrt( ephem(:,2).*ephem(:,2) + ephem(:,3).*ephem(:,3) + ephem(:,4).*ephem(:,4) );
   Rbar_nom   = mean(R)
   Rsigma_nom = std(R)
   diff       =  [Rbar_nom - radius ; Rsigma_nom - 0 ]   

   counter = counter + 1;
end

dos(['copy ephem.txt ',filename]);