function state = state_from_elem(elem, mu)

% Extract elements:
a    = elem(1);  
e    = elem(2);
incl = elem(3);
RAAN = elem(4);
AOP  = elem(5);
M    = elem(6);

% Compute true anomaly:
E = kepler_E(e,M);

    % Employ geometric relationships:
    cos_nu = (cos(E) - e)/(1 - e*cos(E));
    sin_nu = sin(E)*((1 - e^2)^0.5)/(1 - e*cos(E));

    % Run ATAN2 to eliminate quadrant ambiguity:
    nu = atan2(sin_nu,cos_nu);
    
    % [-pi,pi] -> [0,2*pi]
    if nu < 0
        nu = nu + 2*pi;
    end

% Compute state vectors within perifocal frame:
r_p = (a*(1 - e^2))/(1 + e*cos(nu))*[cos(nu); sin(nu); 0];
v_p = sqrt(mu/(a*(1 - e^2)))*[-sin(nu); (e + cos(nu)); 0];

% Compute transformation matrices:
R1 = [ cos(RAAN) sin(RAAN) 0
      -sin(RAAN) cos(RAAN) 0
       0         0         1];

R2 = [ 1  0         0
       0  cos(incl) sin(incl)
       0 -sin(incl) cos(incl)];

R3 = [ cos(AOP) sin(AOP) 0
      -sin(AOP) cos(AOP) 0
       0        0        1];

% Transform frames (Perifocal -> Geocentric ECI):
Q_pX = (R3*R2*R1).';

% Transform state vectors:
r = (Q_pX*r_p);
v = (Q_pX*v_p);

% Output state vectors:
state = [r; v];

end
