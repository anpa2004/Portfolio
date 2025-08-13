function var_dot = AircraftEOMBRYCEHELP(time, aircraft_state, aircraft_surfaces, wind_inertial, aircraft_parameters)
% Computes state derivatives

% Devector
x = aircraft_state(1);         % m
y = aircraft_state(2);         % m
z = aircraft_state(3);         % m
phi = aircraft_state(4);       % rad
theta = aircraft_state(5);     % rad
psi = aircraft_state(6);       % rad
u = aircraft_state(7);         % m/s
v = aircraft_state(8);         % m/s
w = aircraft_state(9);         % m/s
p = aircraft_state(10);        % rad/s
q = aircraft_state(11);        % rad/s
r = aircraft_state(12);        % rad/s

% Calculate Density
density = stdatmo(-z);

% Aircraft Forces and Moments
[aero_forces, aero_moments] = AeroForcesAndMoments(aircraft_state, aircraft_surfaces, wind_inertial, density, aircraft_parameters);

% Aerodynamic Forces
X = aero_forces(1);       % N
Y = aero_forces(2);       % N
Z = aero_forces(3);       % N

% Rotation Moments
L = aero_moments(1);        % N*m
M = aero_moments(2);        % N*m
N = aero_moments(3);        % N*m

% Moments of Inertia
Ix = aircraft_parameters.Ix;
Iy = aircraft_parameters.Iy;
Iz = aircraft_parameters.Iz;
Ixz = aircraft_parameters.Ixz;

% Set Givens
g = aircraft_parameters.g;      % m/s^2
m = aircraft_parameters.m;      % kg

% Calculate Gammas
Gamma = Ix * Iz - Ixz^2;
Gamma1 = (Ixz * (Ix - Iy + Iz)) / Gamma;
Gamma2 = (Iz * (Iz - Iy) + Ixz^2) / Gamma;
Gamma3 = Iz / Gamma;
Gamma4 = Ixz / Gamma;
Gamma5 = (Iz - Ix) / Iy;
Gamma6 = Ixz / Iy;
Gamma7 = (Ix * (Ix - Iy) + Ixz^2) / Gamma;
Gamma8 = Ix / Gamma;

% Derivatives
x_dot = u * cos(theta) * cos(psi) + ...
    v * (sin(phi) * sin(theta) * cos(psi) - cos(phi) * sin(psi)) + ...
    w * (cos(phi) * sin(theta) * cos(psi) + sin(phi) * sin(psi));                   % m/s
y_dot = u * cos(theta) * sin(psi) + ...
    v * (sin(phi) * sin(theta) * sin(psi) + cos(phi) * cos(psi)) + ...
    w * (cos(phi) * sin(theta) * sin(psi) - sin(phi) * cos(psi));                   % m/s
z_dot = u * -sin(theta) + ... 
    v * sin(phi) * cos(theta) + ...
    w * cos(phi) * cos(theta);                                                      % m/s
phi_dot = p * 1 + ...
    q * sin(phi) * tan(theta) + ...
    r * cos(phi) * tan(theta);                                                      % rad/s
theta_dot = p * 0 + ...
    q * cosd(phi) + ...
    r * -sind(phi);                                                                 % rad/s
psi_dot = p * 0 + ...
    q * sin(phi) * sec(theta) + ...
    r * cos(phi) * sec(theta);                                                      % rad/s
u_dot = r * v - q * w + ...
    g * -sin(theta) + ...
    X / m;                                                                          % m/s^2
v_dot = p * w - r * u + ...
    g * cos(theta) * sin(phi) + ...
    Y / m;                                                                          % m/s^2
w_dot = q * u - p * v + ...
    g * cos(theta) * cos(phi) + ...
    Z / m;                                                                          % m/s^2
p_dot = Gamma1 * p * q - Gamma2 * q * r + ...
    Gamma3 * L + Gamma4 * N;                                                        % rad/s^2
q_dot = Gamma5 * p * r - Gamma6 * (p^2 - r^2) + ...
    M / Iy;                                                                         % rad/s^2
r_dot = Gamma7 * p * q - Gamma1 * q * r + ...
    Gamma4 * L + Gamma8 * N;                                                        % rad/s^2

var_dot = [x_dot; y_dot; z_dot; phi_dot; theta_dot; psi_dot; u_dot; v_dot; w_dot; p_dot; q_dot; r_dot];     % [m/s; m/s; m/s; rad/s; rad/s; rad/s; m/s^2; m/s^2; m/s^2; rad/s^2; rad/s^2; rad/s^2]
end