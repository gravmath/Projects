%fp_targeting
function fp_targeter(std_radius,metric,kernel,h,num_pts,slicing,max_iters,filename,derivs_choice)

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
   u = vel(1)
   v = vel(2)
   write_FP_file('\School\Programming\subscribe_only\subscribe_geodesics_parms.txt',M,radius,metric,kernel_choice,h,time,vel,grid,derivs_choice);
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
   Rbar_nom   = Rbar_nom/radius;
   if( Rsigma_nom < best_Rsigma )
      best_u      = u;
      best_v      = v;
      best_Rbar   = Rbar_nom;
      best_Rsigma = Rsigma_nom;
   end
   
   %run du trial
   vel(1)     = u + du;
   vel(2)     = v;
   write_FP_file('\School\Programming\subscribe_only\subscribe_geodesics_parms.txt',M,radius,metric,kernel_choice,h,time,vel,grid,derivs_choice);
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
   ephem      = load('\school\programming\subscribe_only\ephem.txt');
   clear R;
   R          = sqrt( ephem(:,2).*ephem(:,2) + ephem(:,3).*ephem(:,3) + ephem(:,4).*ephem(:,4) );
   Rbar_du    = mean(R)
   Rsigma_du  = std(R)
   diff_du    = [Rbar_du - radius ; Rsigma_du - 0 ]
   Rbar_du    = Rbar_du/radius;
   if( Rsigma_du < best_Rsigma )
      best_u      = vel(1);
      best_v      = vel(2);
      best_Rbar   = Rbar_du;
      best_Rsigma = Rsigma_du;
   end
   
   %run dv trial
   vel(1)     = u;
   vel(2)     = v + dv;
   write_FP_file('\School\Programming\subscribe_only\subscribe_geodesics_parms.txt',M,radius,metric,kernel_choice,h,time,vel,grid,derivs_choice);
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
   Rbar_dv    = mean(R)
   Rsigma_dv  = std(R)   
   diff_dv    = [Rbar_dv - radius ; Rsigma_dv - 0 ]
   Rbar_dv    = Rbar_dv/radius;
   if( Rsigma_dv < best_Rsigma )
      best_u      = vel(1);
      best_v      = vel(2);
      best_Rbar   = Rbar_dv;
      best_Rsigma = Rsigma_dv;
   end
   
   M          = [ (Rbar_du - Rbar_nom)/du    , (Rbar_dv - Rbar_nom)/dv; ...
                  (Rsigma_du - Rsigma_nom)/du, (Rsigma_dv - Rsigma_nom)/dv]
   diff_sc    =  [Rbar_nom - 1; Rsigma_nom - 0 ];
   corrs      = -(M^-1)*diff_sc
   fprintf(fid,'%d\t%20.18f\t%20.18f\t%20.18f\t%20.18f\t%20.18f\t%20.18f\t%20.18f\t%20.18f\t%20.18f\t%20.18f\t%20.18f\t%20.18f\n',counter,u,v,du,dv,Rbar_nom*radius,Rsigma_nom,Rbar_du*radius,Rsigma_du,Rbar_dv*radius,Rsigma_dv,corrs(1),corrs(2));
   vel = [u;v] + corrs
   if ( abs(corrs(1)) < du )
      du = du/5
   end
   if ( abs(corrs(2)) < dv )
      dv = dv/5
   end
   counter = counter + 1;
end

vel(1) = best_u;
vel(2) = best_v;
write_FP_file('\School\Programming\subscribe_only\subscribe_geodesics_parms.txt',M,radius,metric,kernel_choice,h,time,vel,grid,derivs_choice);
%run final nominal
if strcmp(slicing,'bare')
   !subscribe_geodesics_bare
elseif strcmp(slicing,'r3g')
   !subscribe_geodesics_r3g
elseif strcmp(slicing,'spherical')
   !subscribe_geodesics_spherical
elseif strcmp(slicing,'scalar')
   !subscribe_geodesics_scalar
end   

dos(['copy ephem.txt ',filename]);