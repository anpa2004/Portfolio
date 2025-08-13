function xdot = AircraftEOM(time, aircraft_state, aircraft_surfaces,wind_inertial, aircraft_parameters)

var = aircraft_state;

% Devectorizing
phi = var(4);
theta = var(5);
psi = var(6);

uE = var(7);
vE = var(8);
wE = var(9);

p = var(10);
q = var(11);
r = var(12);

% Scalar Inertias from I matrix
Ix = aircraft_parameters.Ix;
Iy = aircraft_parameters.Iy;
Iz = aircraft_parameters.Iz;
Ixz = aircraft_parameters.Ixz;

density = stdatmo(-(aircraft_state(3)));
g = aircraft_parameters.g;
m = aircraft_parameters.m;


v2 = TransformFromBodyToInertial(aircraft_state(7:9),aircraft_state(4:6));
xdot(1)=v2(1);
xdot(2) = v2(2);
xdot(3) = v2(3);

% Taking aero forces/ moments from external functions
[aero_forces, aero_moments] = AeroForcesAndMoments(aircraft_state,aircraft_surfaces, wind_inertial, density, aircraft_parameters);
L = aero_moments(1);
M = aero_moments(2);
N = aero_moments(3);

X = aero_forces(1);
Y = aero_forces(2);
Z = aero_forces(3);

% Gammas
Gamma = (Ix*Iz-Ixz^2);
Gamma1 = Ixz*(Ix-Iy+Iz)/Gamma;
Gamma2 = (Iz*(Iz-Iy)+Ixz^2)/Gamma;
Gamma3 = Iz/Gamma;
Gamma4 = Ixz/Gamma;
Gamma5 = (Iz-Ix)/Iy;
Gamma6 = Ixz/Iy;
Gamma7 = (Ix*(Ix-Iy)+Ixz^2)/Gamma;
Gamma8 = Ix/Gamma;


% Vector EOM 2
A2 = [1 sin(phi)*tan(theta) cos(phi)*tan(theta);0 cos(phi) -sin(phi); 0 sin(phi)*sec(theta) cos(phi)*sec(theta)];
v3 = [p;q;r];
v4 = A2*v3;

xdot(4) = v4(1);
xdot(5) = v4(2);
xdot(6) = v4(3);

% Vector EOM 3
xdot(7) = r*vE-q*wE -g*sin(theta) + X/m;
xdot(8) = p*wE - r*uE + g*cos(theta)*sin(phi) + Y/m;
xdot(9) = q*uE - p*vE +g*cos(theta)*cos(phi) + Z/m;

% Vector EOM 4
xdot(10) = Gamma1*p*q - Gamma2*q*r + Gamma3*L + Gamma4*N;
xdot(11) = Gamma5*p*r - Gamma6*(p^2-r^2) + M/Iy;
xdot(12) = Gamma7*p*q - Gamma1*q*r + Gamma4*L + Gamma8*N;

xdot = xdot';