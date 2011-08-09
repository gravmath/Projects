function E = kepler_E(e,M)

%-------------------------------------------------------------------------
% Programmer: Howard Curtis
%
% This function uses Newton’s method to solve Kepler’s
% equation E - e*sin(E) = M for the eccentric anomaly,
% given the eccentricity and the mean anomaly.
%
% E - eccentric anomaly (radians)
% e - eccentricity, passed from the calling program
% M - mean anomaly (radians), passed from the calling program
% ------------------------------------------------------------------------

% Set error tolerance:
error = 1e-9;

% Set initial value of E:
if M < pi
    E = M + e/2;
else
    E = M - e/2;
end

% Iterate until E is determined within tolerance:
ratio = 1;

while abs(ratio) > error
    ratio = (E - e*sin(E) - M)/(1 - e*cos(E));
    E = E - ratio;
end

end

%-------------------------------------------------------------------------