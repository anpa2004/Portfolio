function PlotAircraftSim(time, aircraft_state_array, control_input_array, fig, col)

% Assign the state array to vectors
x_position = aircraft_state_array(:,1);
y_position = aircraft_state_array(:,2);
z_position = aircraft_state_array(:,3);
roll_angle = rad2deg(aircraft_state_array(:,4));
pitch_angle = rad2deg(aircraft_state_array(:,5));
yaw_angle = rad2deg(aircraft_state_array(:,6));
x_velocity = aircraft_state_array(:,7);
y_velocity = aircraft_state_array(:,8);
z_velocity = aircraft_state_array(:,9);
roll_rate = rad2deg(aircraft_state_array(:,10));
pitch_rate = rad2deg(aircraft_state_array(:,11));
yaw_rate = rad2deg(aircraft_state_array(:,12));

% Expand the control input array
delta_e = linspace(control_input_array(1),control_input_array(1),length(time));
delta_e = rad2deg(delta_e);
delta_a = linspace(control_input_array(2),control_input_array(2),length(time));
delta_a = rad2deg(delta_a);
delta_r = linspace(control_input_array(3),control_input_array(3),length(time));
delta_r = rad2deg(delta_r);
delta_t = linspace(control_input_array(4),control_input_array(4),length(time));

% Create figure 1
figure(fig(1));
sgtitle('Inertial Positions')
subplot(3,1,1)
hold on;
plot(time,x_position, col)
xlabel('Time [s]')
ylabel('X [m]')
hold off;
subplot(3,1,2)
hold on;
plot(time,y_position, col)
xlabel('Time [s]')
ylabel('Y [m]')
hold off;
subplot(3,1,3)
hold on;
plot(time,z_position, col)
xlabel('Time [s]')
ylabel('Z [m]')
hold off;

% Create figure 2
figure(fig(2));
sgtitle('Rotation Angles')
subplot(3,1,1)
hold on;
plot(time,roll_angle, col)
xlabel('Time [s]')
ylabel('Roll Angle [deg]')
hold off;
subplot(3,1,2)
hold on;
plot(time,pitch_angle, col)
xlabel('Time [s]')
ylabel('Pitch Angle [deg]')
hold off;
subplot(3,1,3)
hold on;
plot(time,yaw_angle, col)
xlabel('Time [s]')
ylabel('Yaw Angle [deg]')
hold off;

% Create figure 3
figure(fig(3));
sgtitle('Inertial Velocities')
subplot(3,1,1)
hold on;
plot(time,x_velocity, col)
xlabel('Time [s]')
ylabel('uE [m/s]')
hold off;
subplot(3,1,2)
hold on;
plot(time,y_velocity, col)
xlabel('Time [s]')
ylabel('vE [m/s]')
hold off;
subplot(3,1,3)
hold on;
plot(time,z_velocity, col)
xlabel('Time [s]')
ylabel('wE [m/s]')
hold off;

% Create figure 4
figure(fig(4));
sgtitle('Rotation Rates')
subplot(3,1,1)
hold on;
plot(time,roll_rate, col)
xlabel('Time [s]')
ylabel('p [deg/s]')
hold off;
subplot(3,1,2)
hold on;
plot(time,pitch_rate, col)
xlabel('Time [s]')
ylabel('q [deg/s]')
hold off;
subplot(3,1,3)
hold on;
plot(time,yaw_rate, col)
xlabel('Time [s]')
ylabel('r [deg/s]')
hold off;

% Create figure 5
figure(fig(5))
subplot(4,1,1)
hold on;
plot(time, delta_e', col)
xlabel('Time [s]')
ylabel('\delta_e [N]')
hold off;
subplot(4,1,2)
hold on;
plot(time, delta_a', col)
xlabel('Time [s]')
ylabel('\delta_a [N]')
hold off;
subplot(4,1,3)
hold on;
plot(time, delta_r', col)
xlabel('Time [s]')
ylabel('\delta_r [N]')
hold off;
subplot(4,1,4)
hold on;
plot(time, delta_t', col)
xlabel('Time [s]')
ylabel('\delta_t [N]')
hold off;

% Index the first part of the data values for coloring in figure 6
start = 1:round(3*length(time)/10);
finish = round(9*length(time)/10):round(length(time));
    
% Create figure 6
figure(fig(6))
hold on;
% Done 3 times to account for the line color change
plot3(x_position(:),y_position(:),z_position(:),col)
plot3(x_position(start),y_position(start),z_position(start),'Color','g')
plot3(x_position(finish),y_position(finish),z_position(finish),'Color','r')
title('Aircraft Trajectory')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
legend('Middle','Start','End','Location','southwest')
view(-37.5,30)
hold off;

% Create figure 7
figure(fig(7))
sgtitle('Quadrotor State Values and Control Surfaces')
subplot(8,2,1)
hold on;
plot(time,x_position, col)
xlabel('Time [s]')
ylabel('X [m]')
hold off;
subplot(8,2,3)
hold on;
plot(time,y_position, col)
xlabel('Time [s]')
ylabel('Y [m]')
hold off;
subplot(8,2,5)
hold on;
plot(time,z_position, col)
xlabel('Time [s]')
ylabel('Z [m]')
hold off;
subplot(8,2,7)
hold on;
plot(time,roll_angle, col)
xlabel('Time [s]')
ylabel('Roll [deg]')
hold off;
subplot(8,2,9)
hold on;
plot(time,pitch_angle, col)
xlabel('Time [s]')
ylabel('Pitch [deg]')
hold off;
subplot(8,2,11)
hold on;
plot(time,yaw_angle, col)
xlabel('Time [s]')
ylabel('Yaw [deg]')
hold off;
subplot(8,2,13)
hold on;
plot(time, delta_e', col)
xlabel('Time [s]')
ylabel('\delta_e [deg]')
hold off;
subplot(8,2,15)
hold on;
plot(time, delta_a', col)
xlabel('Time [s]')
ylabel('\delta_a [deg]')
hold off;
subplot(8,2,2)
hold on;
plot(time,x_velocity, col)
xlabel('Time [s]')
ylabel('uE [m/s]')
hold off;
subplot(8,2,4)
hold on;
plot(time,y_velocity, col)
xlabel('Time [s]')
ylabel('vE [m/s]')
hold off;
subplot(8,2,6)
hold on;
plot(time,z_velocity, col)
xlabel('Time [s]')
ylabel('wE [m/s]')
hold off;
subplot(8,2,8)
hold on;
plot(time,roll_rate, col)
xlabel('Time [s]')
ylabel('p [deg/s]')
hold off;
subplot(8,2,10)
hold on;
plot(time,pitch_rate, col)
xlabel('Time [s]')
ylabel('q [deg/s]')
hold off;
subplot(8,2,12)
hold on;
plot(time,yaw_rate, col)
xlabel('Time [s]')
ylabel('r [deg/s]')
hold off;
subplot(8,2,14)
hold on;
plot(time, delta_r', col)
xlabel('Time [s]')
ylabel('\delta_r [deg]')
hold off;
subplot(8,2,16)
hold on;
plot(time, delta_t', col)
xlabel('Time [s]')
ylabel('\delta_t [N/A]')
hold off;
end