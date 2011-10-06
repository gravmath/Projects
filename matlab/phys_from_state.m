function phys = phys_from_state(S, mu)

%unpack the state
R = S(1:3);
V = S(4:6);

% Compute state magnitudes:
r = norm(R);
v = norm(V);

% Compute radial velocity component:
v_r = dot(R,V)/r;

% Compute angular momentum vector/magnitude:
H = cross(R,V);
h = norm(H);

% Compute eccentricity vector/magnitude:
Evec = 1/mu*((v^2 - mu/r)*R - r*(v_r)*V);
e    = norm(Evec);

% Compute SMA:
a = h^2/mu/(1 - e^2);

% Compute the energy
E = -mu/2/a;

elem = [E; H(1); H(2); H(3); Evec(1); Evec(2); Evec(3)];

end
