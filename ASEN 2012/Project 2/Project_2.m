close all; clear; clc;

%Constant Structure
const = getConst();

%time span
tspan = [0,5];

%State Vector
% dV = @(t,y)dVdt_fun(t,y,const);
% [tVair, Vvair] = ode45(dV,tspan,const.v_b);
% figure()
% plot(tVair,Vvair);
% %state = [v_x,a_x,v_z,a_z,-m_w_dot,v_airdt,m_air];

%% Project Verification

verification = load('project2verification.mat');
time_verification = verification.verification.time;
height_verification = verification.verification.height;
distance_verification = verification.verification.distance;
thrust_verification = verification.verification.thrust;
velocity_x_verification = verification.verification.velocity_x;
velocity_y_verification = verification.verification.velocity_y;
volume_air_verification = verification.verification.volume_air;


figure();
plot(distance_verification,height_verification,'LineWidth',1.5);
grid minor;
ylabel('Height (m)');
xlabel('Horizontal distance (m)');
title('Horizontal/Vertical Distance');
ylim([0,max(height_verification)]);

figure();
plot(time_verification,thrust_verification,'LineWidth',1.5);
ylabel('Thrust (N)');
xlabel('Time (s)');
grid minor;
xlim([0,0.5]);
title('Thrust vs Time');

figure();
plot(time_verification,volume_air_verification,'linewidth',1.5);
grid minor;
ylabel('Air Volume (m^3)');
xlabel('Time (s)');
xlim([0,0.5]);

figure();
plot(time_verification,velocity_x_verification,'LineWidth',1.5);
grid minor;
xlabel('Time (s)');
ylabel('X-Velocity (m/s)');
title('X Velocity/ time');

figure();
plot(time_verification,velocity_y_verification,'LineWidth',1.5);
grid minor;
xlabel('Time (s)');
ylabel('Z-Velocity (m/s)');
title('Z Velocity/ time');


%% PASTED 
figure(1);

hold on;
plot(tout,stateOut(:,6));
plot(time_verification,volume_air_verification);
xlim([0,0.5]);

figure(2);
hold on;
plot(tout,thrust);
plot(time_verification,thrust_verification);
xlim([0,.5]);

figure(3);
hold on;
plot(stateOut(:,1),stateOut(:,3),'linewidth',1);
xlabel('x Position (m)');
ylabel('z Position (m)');
grid on;










