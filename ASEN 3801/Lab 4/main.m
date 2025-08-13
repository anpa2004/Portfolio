close all; clear; clc;

var0 = [1;1;-5;0;0;0;0;0;0;0;0;0];
tspan = [0 10];

m = 0.068;
d = 0.06;
km = 0.0024;
Ix = 5.8*10^-5;
Iy = 7.2*10^-5;
Iz = 1*10^-4;
I = [Ix 0 0; 0 Iy 0; 0 0 Iz];
nu = 1*10^-3;
mu = 2*10^-6;

g = 9.81;

force = m*g/4;
motor_forces = [force,force,force,force];

f = @(t,var)QuadrotorEOM(t,var,g,m,I,d,km,nu,mu,motor_forces);

[tvec,var] = ode45(f,tspan,var0);

motor_control = motor_forces(1)*ones(length(tvec),4);
col = 'b-';
fig1 = 1:1:7;
% PlotAircraftSim(tvec,var,motor_control,fig1,col)

%% Lab task 1: 3a

%var0 = [0;0;0;0.03746;0;0;0;4.99649;-0.18726;0;0;0];
w = m*g;
Va = 5;
f_d = nu*Va^2;
phi = atan(f_d/w);
Zc = f_d/sin(phi);
vE = Va*cos(phi);
wE = Va*sin(phi);

var0 = [0;0;0;phi;0;0;0;vE;-wE;0;0;0];

force3 = Zc/4;
motor_forces = [force3,force3,force3,force3];

f3 = @(t,var)QuadrotorEOM(t,var,g,m,I,d,km,nu,mu,motor_forces);

[tvec3,var3] = ode45(f3,tspan,var0);
fig3 = 8:14;
motor_control = motor_forces(1)*ones(length(tvec3),4);
%PlotAircraftSim(tvec3,var3,motor_control,fig3,col)

%% Lab task 1: 3b

Va = 5;
f_d = nu*Va^2;
theta = atan(f_d/w);
Zc = f_d/sin(theta);
uE = Va*cos(theta);
wE = Va*sin(theta);
psi = pi/2;


var0 = [0;0;0;0;-theta;psi;uE;0;-wE;0;0;0];

force3b = Zc/4;
motor_forces = [force3b,force3b,force3b,force3b];

f3b = @(t,var)QuadrotorEOM(t,var,g,m,I,d,km,nu,mu,motor_forces);

[tvec3b,var3b] = ode45(f3b,tspan,var0);
fig3b = 15:21;
motor_control = motor_forces(1)*ones(length(tvec3b),4);
PlotAircraftSim(tvec3b,var3b,motor_control,fig3b,col)

