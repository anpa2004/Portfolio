close all; clear; clc;

aircraft_parameters = ttwistor();
aircraft_state = [0;0;-1609.34;0;0;0;21;0;0;0;0;0];
de = 0;
da = 0;
dr = 0;
dt = 0;
aircraft_surfaces = [de da dr dt]';
wind_inertial = [0;0;0];

density = stdatmo(abs(aircraft_state(3)));

[aero_forces, aero_moments] = AeroForcesAndMoments(aircraft_state,aircraft_surfaces, wind_inertial, density, aircraft_parameters);
aee = aero_forces/aircraft_parameters.m;