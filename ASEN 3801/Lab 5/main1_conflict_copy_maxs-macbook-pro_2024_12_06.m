close all; clear; clc;
col = "b-";
% e,a,r,t
%% Problem 2.1

%Initial State Vector
x0 = [0;0;-1609.34;0;0;0;21;0;0;0;0;0];

%Calling other constants
aircraft_parameters = ttwistor();
de = 0;
da = 0;
dr = 0;
dt = 0;
aircraft_surfaces = [de da dr dt]';
wind_inertial = [0;0;0];
density = stdatmo(abs(x0(3)));

%Finding the motion of the Aircraft
f = @(t,state)AircraftEOM(t,state,aircraft_surfaces,wind_inertial,aircraft_parameters);

tspan = [0 10];
[time,state] = ode45(f,tspan,x0);

%creating the control inputs as arrays for graphing
control_inputs1(:,1) = ones(length(time),1)*de;
control_inputs1(:,2) = ones(length(time),1)*da;
control_inputs1(:,3) = ones(length(time),1)*dr;
control_inputs1(:,4) = ones(length(time),1)*dt;

%Plotting the resulting motion
% figure()
% hold on
% plot3(state(:,1),state(:,2),abs(state(:,3)),'Linewidth',1)
% scatter3(state(1,1),state(1,2),abs(state(1,3)),'g')
% scatter3(state(end,1),state(end,2),abs(state(end,3)),'r')
% grid minor
% xlabel("x (meters)")
% ylabel("y (meters)")
% zlabel("z (meters)")
% view([-29 35])
% hold off
fig1 = 1:7;
%PlotAircraftSim(time,state,control_inputs1,fig1,col);

%% Problem 2.2

%Initial State Vector
x0 = [0,0,-1800,0,0.0278,0,20.99,0,0.5837,0,0,0]';
%Control Surface Deflections
aircraft_surfaces = [0.1079,0,0,0.3182]';

%Simulation of EOM
f = @(t,state)AircraftEOM(t,state,aircraft_surfaces,wind_inertial,aircraft_parameters);

tspan = [0 10];
[time,state] = ode45(f,tspan,x0);

%Control Inputs as Arrays for Graphing
control_inputs2(:,1) = ones(length(time),1)*aircraft_surfaces(1);
control_inputs2(:,2) = ones(length(time),1)*aircraft_surfaces(2);
control_inputs2(:,3) = ones(length(time),1)*aircraft_surfaces(3);
control_inputs2(:,4) = ones(length(time),1)*aircraft_surfaces(4);

% figure()
% hold on
% plot3(state(:,1),state(:,2),abs(state(:,3)),'Linewidth',1)
% scatter3(state(1,1),state(1,2),abs(state(1,3)),'g')
% scatter3(state(end,1),state(end,2),abs(state(end,3)),'r')
% grid minor
% xlabel("x (meters)")
% ylabel("y (meters)")
% zlabel("z (meters)")
% view([-29 35])
% hold off
fig2 = 8:14;
%PlotAircraftSim(time,state,control_inputs2,fig2,col);

%% Problem 2.3

%deg to rad conversion
c = pi/180;
%initial state vector
x0 = [0,0,-1800,15*c,-12*c,270*c,19,3,-2,0.08*c,-0.2*c,0]';
%control surfaces
aircraft_surfaces = [5*c;2*c;-13*c;0.3];

%simulating EOM
f = @(t,state)AircraftEOM(t,state,aircraft_surfaces,wind_inertial,aircraft_parameters);
tspan = [0 200]; %[s]
[time,state] = ode45(f,tspan,x0); % ERROR IN INT TOLERANCES

% figure()
% hold on
% plot3(state(:,1),state(:,2),abs(state(:,3)),'Linewidth',1)
% scatter3(state(1,1),state(1,2),abs(state(1,3)),'g')
% scatter3(state(end,1),state(end,2),abs(state(end,3)),'r')
% grid minor
% xlabel("x (meters)")
% ylabel("y (meters)")
% zlabel("z (meters)")
% view([-29 35])
% hold off
fig3 = 15:21;
control_inputs3(:,1) = ones(length(time),1)*aircraft_surfaces(1);
control_inputs3(:,2) = ones(length(time),1)*aircraft_surfaces(2);
control_inputs3(:,3) = ones(length(time),1)*aircraft_surfaces(3);
control_inputs3(:,4) = ones(length(time),1)*aircraft_surfaces(4);
%PlotAircraftSim(time,state,control_inputs3,fig3,col);

%% Problem 3.1
clc; 

c = pi/180;
%Initial Conditions
tspan = [0 100];
x0 = [0,0,-1800,15*c,-12*c,279*c,19,3,-2,0.08*c,-0.2*c,0];
aircraft_surfaces = [0.1079,0,0,0.3182]';

%Doublet 
doublet_size = 25*c;
doublet_time = 1;

f = @(time,aircraft_state)AircraftEOMDoublet(time, aircraft_state, aircraft_surfaces,doublet_size, doublet_time, wind_inertial, aircraft_parameters);

[time,state] = ode45(f,tspan,x0);

for i = 1:length(time)
    control_surfaces(i,:) = AircraftEOMDoubletGETDELTAS(time(i),state(i,:)',aircraft_surfaces,doublet_size,doublet_time,wind_inertial,aircraft_parameters)';
end

% figure()
% hold on
% plot3(state(:,1),state(:,2),abs(state(:,3)),'Linewidth',1)
% scatter3(state(1,1),state(1,2),abs(state(1,3)),'g')
% scatter3(state(end,1),state(end,2),abs(state(end,3)),'r')
% grid minor
% xlabel("x (meters)")
% ylabel("y (meters)")
% zlabel("z (meters)")
% view([-29 35])
% hold off
%PlotAircraftSim(time,state,)

fig4 = 22:28;

PlotAircraftSim(time,state,control_surfaces,fig4,col)