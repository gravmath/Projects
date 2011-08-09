function elem = elem_from_state(R, V, mu)

% Compute state magnitudes:
r = norm(R);
v = norm(V);

% Compute radial velocity component:
v_r = dot(R,V)/r;

% Compute angular momentum vector/magnitude:
H = cross(R,V);
h = norm(H);

% Compute inclination:
incl = acos(H(3)/h);

% Compute node vector/magnitude:
N = cross([0 0 1],H);
n = norm(N);

% Compute RAAN (quandrant ambiguity):
if n ~= 0
    RAAN = acos(N(1)/n);
    if N(2) < 0
        RAAN = 2*pi - RAAN;
    end
else
    RAAN = 0;
end

% Compute eccentricity vector/magnitude:
E = 1/mu*((v^2 - mu/r)*R - r*(v_r)*V);
e = norm(E);

% Compute AOP (quandrant ambiguity):
if n ~= 0
    if e > 1e-9
        AOP = acos(dot(N,E)/n/e);
        if E(3) < 0
           AOP = 2*pi - AOP;
        end
    else
        AOP = 0;
    end
else
    AOP = 0;
end

% Compute true anomaly (quadrant ambiguity):
if e > 1e-9
    nu = acos(dot(E,R)/e/r);
    if v_r < 0
        nu = 2*pi - nu;
    end
else
    cp = cross(N,R);
    if cp(3) >= 0
        nu = acos(dot(N,R)/n/r);
    else
        nu = 2*pi - acos(dot(N,R)/n/r);
    end
end

% Compute SMA:
a = h^2/mu/(1 - e^2);

% Compute mean anomaly:
cos_E = (e + cos(nu))/(1 + e*cos(nu));
sin_E = (sin(nu)*sqrt(1 - e^2))/(1 + e*cos(nu));
E = atan2(sin_E,cos_E);
M = E - e*sin(E);

% Output Keplerian elements:
elem = [a; e; incl; RAAN; AOP; M];

end
