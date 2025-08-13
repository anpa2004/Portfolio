% Housekeeping
close;
clear;
clc;

% Givens
wind = [0; 0; 0];

% Call Aircraft Parameters
aircraft_parameters = ttwistor;

%% Problem 2
col = 'b-';


%% Problem 2.3
% Set State Vectors

x2c = [0; 0; -1800; deg2rad(15); deg2rad(-12); deg2rad(270); 19; 3; -2; deg2rad(0.08); deg2rad(-0.2); 0];
u2c = [deg2rad(5); deg2rad(2); deg2rad(-13); 0.3];

% Set Time Span
tspan = [0 200];

% Call EOM
[time2c, state2c] = ode45(@(t,x2c)AircraftEOMBRYCEHELP(t,x2c,u2c,wind,aircraft_parameters),tspan,x2c);

% Plot
fig = 15:21;
PlotAircraftSimBRYCEHELP(time2c, state2c, u2c, fig, col);

%% Problem 3
% Doublet Parameters
doublet_size = deg2rad(25);
doublet_time = 1;

%% Problem 3.1
% Set Time Span
tspan = [0, 4];

% EOM
[time3a, state3a] = ode45(@(t,x2b)AircraftEOMBRYCEHELPDoublet(t,x2b,u2b,doublet_size,doublet_time,wind,aircraft_parameters),tspan,x2b);

% Plot
fig = 22:28;
PlotAircraftSim(time3a, state3a, u2b, fig, col);

%% Problem 3.2
% Set Time Span
tspan = [0, 100];

% EOM
[time3b, state3b] = ode45(@(t,x2b)AircraftEOMBRYCEHELPDoublet(t,x2b,u2b,doublet_size,doublet_time,wind,aircraft_parameters),tspan,x2b);

% Plot
fig = 29:35;
PlotAircraftSim(time3b, state3b, u2b, fig, col);