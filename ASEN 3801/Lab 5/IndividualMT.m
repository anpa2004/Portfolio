clc;
close all;
clear all;

%Constants
h = 1609.34;    %m
Va = 21;    %m/s
rho = stdatmo(h);

wind_E = [0;0;0];

aircraft_parameters = ttwistor();

%State vector
X = [0;0;-h;0;0;0;Va;0;0;0;0;0];
controls = [0;0;0;0];

[aero_forces,aero_moments] = AeroForcesAndMoments(X, controls, wind_E, rho, aircraft_parameters);

%Calculating inertial acceleration in inertial coordinates
aE = aero_forces ./ aircraft_parameters.m + [0;0;aircraft_parameters.g]

