%-------------------------------------------------------------------------
% NUMERICAL INTEGRATION OF THE TWO-BODY PROBLEM WITH J2 PERTURBATION
% Programmer: Paul V. Anderson, Navigation and Mission Design Branch
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
clc, clear all, close all;
delete('kepler_time.txt');
delete('kepler_state.txt');
tic;
%-------------------------------------------------------------------------

% Define necessitated numerical constants:
mu         = 398600.4418;  % [km^3/sec^2]
RE         = 6378.137;     % [km]
J2         = 1.0826269E-3; % Vallado page 521 - possible JGM-3
J2         = 1.0826400E-3; % BMW page 422
num_orbits = 30; 

          % [a; e; incl; RAAN; AOP; M]
init_cond = [42095; 0.8181818; 28.5*(pi/180); 0; 0; 0];

          % [Rx; Ry; Rz; Vx; Vy; Vz]
init_cond = state_from_elem(init_cond, mu);

% Define simulation parameters:
r         = init_cond(1:3);
v         = init_cond(4:6);
E         = 1/2*(v'*v) - mu/sqrt(r'*r);
SMA       = -mu/2/E;
P         = 2*pi*sqrt(SMA^3/mu);        
time_span = [0, 600];

% Numerical integration of KEPLER.m model:
options = odeset('AbsTol', 1e-11, 'RelTol', 1e-11);
[t,S]   = ode45('KEPLER', time_span, init_cond, options, mu, J2, RE);

% Output/print results to files:
fid = fopen('kepler_time.txt','a');
for n = 1:length(t)
    fprintf(fid, '%E\n', t(n));
end
fclose(fid); 

fid = fopen('kepler_state.txt','a');
for m = 1:length(t)
    for n = 1:6
        fprintf(fid, '%E\t', S(m,n));
    end
    fprintf(fid,'\n');
end
fclose(fid);

% Extract position/velocity data:
S_pos = S(:,1:3);
S_vel = S(:,4:6);

% Compute position/velocity magnitudes:
r_mag = sqrt(sum(S_pos.*S_pos,2));
v_mag = sqrt(sum(S_vel.*S_vel,2));

% Graphical representations:
figure(1)
subplot(211)
plot(t, r_mag, '-r', 'LineWidth', 2);
title('Variation in Orbital Radius with Simulation Time','FontWeight','b');
xlabel('Simulation Time (sec)');
ylabel('Orbital Radius (km)');
axis tight
grid on

subplot(212)
plot(t, v_mag, '-r', 'LineWidth', 2);
title('Variation in Orbital Velocity with Simulation Time','FontWeight','b');
xlabel('Simulation Time (sec)');
ylabel('Orbital Velocity (km/sec)');
axis tight
grid on

% Preallocation of elem array:
elem = zeros(length(t),6);

% Plot variation in Keplerian elements:
for m = 1:length(t)
    elem(m,:) = elem_from_state(S_pos(m,:), S_vel(m,:), mu);
end

figure(2)
plot(t, elem(:,1), 'b-', 'LineWidth', 2);
title('Variation in SMA with Simulation Time','FontWeight','b');
xlabel('Simulation Time (sec)');
ylabel('SMA (km)');
grid on

figure(3)
plot(t, elem(:,2), 'b-', 'LineWidth', 2);
title('Variation in Eccentricity with Simulation Time','FontWeight','b');
xlabel('Simulation Time (sec)');
ylabel('Eccentricity');
grid on

figure(4)
plot(t, elem(:,3), 'b-', 'LineWidth', 2);
title('Variation in Inclination with Simulation Time','FontWeight','b');
xlabel('Simulation Time (sec)');
ylabel('Inclination (rad)');
grid on

figure(5)
plot(t, elem(:,4), 'b-', 'LineWidth', 2);
title('Variation in RAAN with Simulation Time','FontWeight','b');
xlabel('Simulation Time (sec)');
ylabel('RAAN (rad)');
grid on

figure(6)
plot(t, elem(:,5), 'b-', 'LineWidth', 2);
title('Variation in AOP with Simulation Time','FontWeight','b');
xlabel('Simulation Time (sec)');
ylabel('AOP (rad)');
grid on

figure(7)
plot(t, elem(:,6), 'b-', 'LineWidth', 2);
title('Variation in Mean Anomaly with Simulation Time','FontWeight','b');
xlabel('Simulation Time (sec)');
ylabel('Mean Anomaly (rad)');
grid on

% Plot planar orbital trajectory:
figure(8)    
plot(S_pos(:,1), S_pos(:,2),'Color',[1,0,0],'LineWidth',2);
title('2D Orbital Trajectory','FontWeight','b');
xlabel('Abscissa (km)');
ylabel('Ordinate (km)');
axis square
grid on

figure(9)

%{
    % Plot spherical Earth representation:
    phi_earth = linspace(0,pi);
    theta_earth = linspace(0,2*pi);
    [phi_earth,theta_earth] = meshgrid(phi_earth,theta_earth);

    x_earth = RE*sin(phi_earth).*cos(theta_earth);
    y_earth = RE*sin(phi_earth).*sin(theta_earth);
    z_earth = RE*cos(phi_earth);
    
    surf(x_earth,y_earth,z_earth);
    colormap winter, shading interp, alpha(0.3)
    set(gca,'DataAspectRatio',[1,1,1],'Color',[0,0,0.2]);
    hold on
%}

    % Plot Earth globe:
    load topo
    [x,y,z] = sphere(50);
    cla reset, axis square

    props.AmbientStrength = 0.1;
    props.DiffuseStrength = 1;
    props.SpecularColorReflectance = 0.5;
    props.SpecularExponent = 20;
    props.SpecularStrength = 1;
    props.FaceColor= 'texture';
    props.EdgeColor = 'none';
    props.FaceLighting = 'phong';
    props.Cdata = topo;

    surface(RE*x, RE*y, RE*z, props);
    view(3), camzoom(1.5)
    set(gca,'DataAspectRatio',[1,1,1],'Color',[0,0,0]);
    hold on

% Plot orbital trajectory:
plot3(S_pos(:,1), S_pos(:,2), S_pos(:,3),'-w','LineWidth',1.5);
title('3D Orbital Trajectory','FontWeight','b');
xlabel('Abscissa (km)');
ylabel('Ordinate (km)');
zlabel('Applicate (km)');
grid on
hold off

%-------------------------------------------------------------------------
toc;
%-------------------------------------------------------------------------

figure(10);

% Animate satellite and orbital trajectory:
for p = 1:(0.2*length(t))
    
%{
    % Plot spherical Earth representation:
    phi_earth = linspace(0,pi);
    theta_earth = linspace(0,2*pi);
    [phi_earth,theta_earth] = meshgrid(phi_earth,theta_earth);

    x_earth = RE*sin(phi_earth).*cos(theta_earth);
    y_earth = RE*sin(phi_earth).*sin(theta_earth);
    z_earth = RE*cos(phi_earth);
    
    surf(x_earth,y_earth,z_earth);
    colormap winter, shading interp, alpha(0.3), camzoom(1.5)
    set(gca,'DataAspectRatio',[1,1,1],'Color',[0,0,0.2]);
    hold on
%}    
    
    % Plot Earth globe:
    load topo
    [x,y,z] = sphere(50);
    cla reset, axis square

    props.AmbientStrength = 0.1;
    props.DiffuseStrength = 1;
    props.SpecularColorReflectance = 0.5;
    props.SpecularExponent = 20;
    props.SpecularStrength = 1;
    props.FaceColor= 'texture';
    props.EdgeColor = 'none';
    props.FaceLighting = 'phong';
    props.Cdata = topo;

    surface(RE*x, RE*y, RE*z, props);
    view(3), camzoom(2)
    set(gca,'DataAspectRatio',[1,1,1],'Color',[0,0,0]);
    set(gcf,'Color',[0,0,0]);
    hold on
    
    % Plot representation of satellite and velocity vector:
    plot3(S_pos(5*p,1), S_pos(5*p,2), S_pos(5*p,3),'-r*','MarkerEdgeColor','r','MarkerSize',8);
    quiver3(S(5*p,1),S(5*p,2),S(5*p,3),S(5*p,4),S(5*p,5),S(5*p,6),200,'-m','LineWidth',1.5);
    
    %{
    
    % Change color of animated trajectory each period completed:
    if t(5*p) >= P
        for m = 1:num_orbits
            rev = find(t >= m*P);
            if t(5*p) >= m*P && t(5*p) < (m + 1)*P
               plot3(S_pos(rev(1):5*p,1), S_pos(rev(1):5*p,2), S_pos(rev(1):5*p,3),'Color',[0.03*m,0.03*m,0.03*m],'LineWidth',1.5);
               title('3D Orbital Trajectory','FontWeight','b');
               xlabel('Abscissa (km)');
               ylabel('Ordinate (km)');
               zlabel('Applicate (km)');
               grid on
               hold on
               plot3(S_pos(1:rev(1),1), S_pos(1:rev(1),2), S_pos(1:rev(1),3),'-w','LineWidth',1.5);
               title('3D Orbital Trajectory','FontWeight','b');
               xlabel('Abscissa (km)');
               ylabel('Ordinate (km)');
               zlabel('Applicate (km)');
               grid on
               hold off
               M(5*p) = getframe;
            end
        end   
    else
        plot3(S_pos(1:5*p,1), S_pos(1:5*p,2), S_pos(1:5*p,3),'-w','LineWidth',1.5);
        title('3D Orbital Trajectory','FontWeight','b');
        xlabel('Abscissa (km)');
        ylabel('Ordinate (km)');
        zlabel('Applicate (km)');
        grid on
        hold off
        M(5*p) = getframe;
    end
    
    %}     
    
    % Plot orbital trajectory:
    plot3(S_pos(1:5*p,1), S_pos(1:5*p,2), S_pos(1:5*p,3),'-w','LineWidth',1.5);
    title('3D Orbital Trajectory','FontWeight','b');
    hold off
    
    % Obtain animation frames:
    M(p) = getframe;
    
end

% Play animation:
movie(M,1,20);

%-------------------------------------------------------------------------
% END OF CODE
%-------------------------------------------------------------------------