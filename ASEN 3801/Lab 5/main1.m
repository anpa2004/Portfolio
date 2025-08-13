%% Main

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

tspan = [0 200];
[time,state] = ode45(f,tspan,x0);

%creating the control inputs as arrays for graphing
control_inputs1(:,1) = ones(length(time),1)*de;
control_inputs1(:,2) = ones(length(time),1)*da;
control_inputs1(:,3) = ones(length(time),1)*dr;
control_inputs1(:,4) = ones(length(time),1)*dt;

fig1 = 1:7;
PlotAircraftSim(time,state,control_inputs1,fig1,col);

%% Problem 2.2

%Initial State Vector
x0 = [0,0,-1800,0,0.0278,0,20.99,0,0.5837,0,0,0]';
%Control Surface Deflections
aircraft_surfaces = [0.1079,0,0,0.3182]';

%Simulation of EOM
f = @(t,state)AircraftEOM(t,state,aircraft_surfaces,wind_inertial,aircraft_parameters);

tspan = [0 200];
[time,state] = ode45(f,tspan,x0);

%Control Inputs as Arrays for Graphing
control_inputs2(:,1) = ones(length(time),1)*aircraft_surfaces(1);
control_inputs2(:,2) = ones(length(time),1)*aircraft_surfaces(2);
control_inputs2(:,3) = ones(length(time),1)*aircraft_surfaces(3);
control_inputs2(:,4) = ones(length(time),1)*aircraft_surfaces(4);

fig2 = 8:14;
PlotAircraftSim(time,state,control_inputs2,fig2,col);

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

fig3 = 15:21;
control_inputs3(:,1) = ones(length(time),1)*aircraft_surfaces(1);
control_inputs3(:,2) = ones(length(time),1)*aircraft_surfaces(2);
control_inputs3(:,3) = ones(length(time),1)*aircraft_surfaces(3);
control_inputs3(:,4) = ones(length(time),1)*aircraft_surfaces(4);
PlotAircraftSim(time,state,control_inputs3,fig3,col);

%% Problem 3.1
clc; 

%Initial Conditions
tspan = [0 4];
x0 = [0,0,-1800,0,0.0278,0,20.99,0,0.5837,0,0,0]';
aircraft_surfaces = [0.1079,0,0,0.3182]';

%Doublet 
doublet_size = 25*c;
doublet_time = 1;

%Simulating EOM
f = @(time,aircraft_state)AircraftEOMDoublet(time, aircraft_state, aircraft_surfaces,doublet_size, doublet_time, wind_inertial, aircraft_parameters);

[time,state] = ode45(f,tspan,x0);

for i = 1:length(time)
    control_surfaces(i,:) = AircraftEOMDoubletGETDELTAS(time(i),state(i,:)',aircraft_surfaces,doublet_size,doublet_time,wind_inertial,aircraft_parameters)';
end

fig4 = 22:28;

PlotAircraftSim(time,state,control_surfaces,fig4,col)

%% Analysis of 3.1 Results
pt_sp = find(time>2); %short period mode points

max_val(1) = max(state(pt_sp,11));
max_val(2) = max(state(124:end,11));

pt_max(1) = find(state(:,11)==max_val(1));
pt_max(2) = find(state(:,11)==max_val(2));


%period 
T = diff(time(pt_max)); %the time between two peaks/troughs
wd = 2*pi/T;

%logmarithic decrement
delta = log(max_val(1)/max_val(2));
zeta = delta/sqrt(4*pi^2+delta^2);

%damping and natural frequency
wn = wd/sqrt(1-zeta^2);

%% Problem 3.2
clc; 

%Initial Conditions
tspan = [0 100];
x0 = [0,0,-1800,0,0.0278,0,20.99,0,0.5837,0,0,0]';
aircraft_surfaces = [0.1079,0,0,0.3182]';

%Doublet 
doublet_size = 25*c;
doublet_time = 1;

%Simulating EOM
f = @(time,aircraft_state)AircraftEOMDoublet(time, aircraft_state, aircraft_surfaces,doublet_size, doublet_time, wind_inertial, aircraft_parameters);

[time,state] = ode45(f,tspan,x0);

[zeta,omegan] = FindDampingNaturalFreq(time,state);
for i = 1:length(time)
    control_surfaces(i,:) = AircraftEOMDoubletGETDELTAS(time(i),state(i,:)',aircraft_surfaces,doublet_size,doublet_time,wind_inertial,aircraft_parameters)';
end

fig4 = 29:35;
ts = 3/(zeta*omegan)
PlotAircraftSim(time,state,control_surfaces,fig4,col)