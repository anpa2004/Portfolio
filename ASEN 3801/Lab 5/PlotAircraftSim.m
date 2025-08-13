function PlotAircraftSim(time, aircraft_state_array, control_input_array, fig, col)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function plots the aircraft state vector vs time, aircraft
% trajectory and the forces vs time on subplots and on one overview figure.
% It also overlays the lineraized and non lineraized dynamics to compare them 
%
% OUTPUTS: 
% None 
% 
% INPUTS: 
% time: time 
% aircraft_state_array = state changing over time
% fig: figure numbers array 
% col: color for plotting  
% isOverlay: condition (true or false) to overlay linerarized dynamics onto
% nonlinear dynamics plots 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if nargin < 6 % Number of aircraft sim variables, so if less then 6 overlay
%     %not required
%     isOverlay = false; % Default is to not overlay.
% end

% Inertial Position

figure(fig(1))
%inertial x position
subplot(3,1,1)
hold on;
plot(time,aircraft_state_array(:,1),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("Inertial x Position (m)");
%inertial y position
subplot(3,1,2)
hold on;
plot(time,aircraft_state_array(:,2),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("Inertial y Position (m)");
%inertial z position
subplot(3,1,3)
hold on;
plot(time,aircraft_state_array(:,3),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("Inertial z Position (m)");
sgtitle("Inertial Position Componentwise Graphs")

% Euler Angles Graph
figure(fig(2))
%phi angle
subplot(3,1,1)
hold on;
plot(time,aircraft_state_array(:,4),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("\phi (rad)");
%theta angle
subplot(3,1,2)
hold on;
plot(time,aircraft_state_array(:,5),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("\theta (rad)");
%psi angle
subplot(3,1,3)
hold on;
plot(time,aircraft_state_array(:,6),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("\psi (rad)");
sgtitle("Euler Angles")

% Inertial Velocity in Body frame
figure(fig(3))
%u-velocity
subplot(3,1,1)
hold on;
plot(time,aircraft_state_array(:,7),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("Inertial x Velocity (m/s)");
%v-velocity
subplot(3,1,2)
hold on;
plot(time,aircraft_state_array(:,8),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("Inertial y Velocity (m/s)");
%w-velocity
subplot(3,1,3)
hold on;
plot(time,aircraft_state_array(:,9),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("Inertial z Velocity (m/s)");

sgtitle("Inertial Velocity in Body Frame")

% Angular velocity
figure(fig(4))
% p-angle rate
subplot(3,1,1)
hold on;
plot(time,aircraft_state_array(:,10),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("Roll Rate (rad/s)");
% q-angle rate
subplot(3,1,2)
hold on;
plot(time,aircraft_state_array(:,11),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("Pitch Rate (rad/s)");
% r-anlge rate
subplot(3,1,3)
hold on;
plot(time,aircraft_state_array(:,12),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("Yaw Rate (rad/s)");

sgtitle("Angular Velocity");

% Control Input Variables
figure(fig(5))
% change in elevation
subplot(4,1,1)
hold on;
plot(time,control_input_array(:,1),col,"Linewidth",1);
grid minor
xlabel("Time (s)");
ylabel("\delta_e")
% change in aileron
subplot(4,1,2)
hold on;
plot(time,control_input_array(:,2),col,"Linewidth",1);
grid minor
xlabel("Time (s)");
ylabel("\delta_a")
% change in rudder
subplot(4,1,3)
hold on;
plot(time,control_input_array(:,3),col,"Linewidth",1);
grid minor
xlabel("Time (s)");
ylabel("\delta_r")
% change in thrust
subplot(4,1,4)
hold on;
plot(time,control_input_array(:,4),col,"Linewidth",1);
grid minor
xlabel("Time (s)");
ylabel("\delta_t")

% Flight Path
figure(fig(6))
hold on;
plot3(aircraft_state_array(:,1),aircraft_state_array(:,2),abs(aircraft_state_array(:,3)),col,'Linewidth',1);
scatter3(aircraft_state_array(1,1),aircraft_state_array(1,2),abs(aircraft_state_array(1,3)),'g')
scatter3(aircraft_state_array(end,1),aircraft_state_array(end,2),abs(aircraft_state_array(end,3)),'r')
grid minor;
xlabel("X Position (m)");
ylabel('Y Position (m)');
zlabel("Z Position (m)");
title("Aircraft Position");
view([-29 35])
hold off;

% Figure 7: All State Components Plus Control Deflections
figure(fig(7))
%x-position
subplot(8,2,1)
hold on;
plot(time,aircraft_state_array(:,1),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("x (m)");
%y-position
subplot(8,2,3)
hold on;
plot(time,aircraft_state_array(:,2),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("y (m)");
%z-position
subplot(8,2,5)
hold on;
plot(time,aircraft_state_array(:,3),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("z (m)");
%phi-angle
subplot(8,2,7)
hold on;
plot(time,aircraft_state_array(:,4),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("\phi (rad)");
%theta-angle
subplot(8,2,9)
hold on;
plot(time,aircraft_state_array(:,5),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("\theta (rad)");
%psi-angle
subplot(8,2,11)
hold on;
plot(time,aircraft_state_array(:,6),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("\psi (rad)");
%
subplot(8,2,2)
hold on;
plot(time,aircraft_state_array(:,7),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("uE (m/s)");
subplot(8,2,4)
hold on;
plot(time,aircraft_state_array(:,8),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("vE (m/s)");
subplot(8,2,6)
hold on;
plot(time,aircraft_state_array(:,9),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("wE (m/s)");
subplot(8,2,8)
hold on;
plot(time,aircraft_state_array(:,10),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("p (rad/s)");
subplot(8,2,10)
hold on;
plot(time,aircraft_state_array(:,11),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("q (rad/s)");
subplot(8,2,12)
hold on;
plot(time,aircraft_state_array(:,12),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("r (rad/s)");
subplot(8,2,13)
hold on;
plot(time,control_input_array(:,1),col,"linewidth",1);
xlabel("Time");
ylabel("\delta_e")
grid minor
subplot(8,2,14)
hold on;
plot(time,control_input_array(:,3),col,"linewidth",1);
xlabel("Time");
ylabel("\delta_r ")
grid minor
subplot(8,2,15)
hold on;
plot(time,control_input_array(:,2),col,"linewidth",1);
xlabel("Time");
ylabel("\delta_a ")
grid minor
subplot(8,2,16)
hold on;
plot(time,control_input_array(:,4),col,"linewidth",1);
xlabel("Time");
ylabel("\delta_t")

grid minor

end