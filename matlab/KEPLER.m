function ds = KEPLER(t, s, dummy, mu, J2, RE)

r = [s(1); s(2); s(3)];
v = [s(4); s(5); s(6)];

% Implement two-body EOM:
r_mag = sqrt(r.'*r);
k = mu/(r_mag^3);

dv_x = -k*s(1)*(1 - J2*(3/2)*((RE/r_mag)^2)*(5*((s(3)/r_mag)^2) - 1));
dv_y = (s(2)/s(1))*dv_x;
dv_z = -k*s(3)*(1 + J2*(3/2)*((RE/r_mag)^2)*(3 - 5*((s(3)/r_mag)^2)));

% Compute differential state vector:
ds = [s(4); s(5); s(6); dv_x; dv_y; dv_z];

end