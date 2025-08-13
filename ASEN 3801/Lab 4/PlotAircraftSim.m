function PlotAircraftSim(time,aircraft_state_array,control_input_array,fig,col)

% Inertial Position
figure(fig(1))
subplot(3,1,1)
plot(time,aircraft_state_array(:,1),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("Inertial x Position (m)");
subplot(3,1,2)
plot(time,aircraft_state_array(:,2),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("Inertial y Position (m)");
subplot(3,1,3)
plot(time,aircraft_state_array(:,3),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("Inertial z Position (m)");
sgtitle("Inertial Position Componentwise Graphs")

% Euler Angles Graph
figure(fig(2))
subplot(3,1,1)
plot(time,aircraft_state_array(:,4),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("Euler Angle \phi (rad)");
subplot(3,1,2)
plot(time,aircraft_state_array(:,5),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("Euler Angle \theta (rad)");
subplot(3,1,3)
plot(time,aircraft_state_array(:,6),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("Euler Angle \psi (rad)");
sgtitle("Euler Angles")

% Inertial Velocity in Body frame
figure(fig(3))
subplot(3,1,1)
plot(time,aircraft_state_array(:,7),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("Inertial x Velocity (m/s)");
subplot(3,1,2)
plot(time,aircraft_state_array(:,8),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("Inertial y Velocity (m/s)");
subplot(3,1,3)
plot(time,aircraft_state_array(:,9),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("Inertial z Velocity (m/s)");

sgtitle("Inertial Velocity in Body Frame")

% Angular velocity
figure(fig(4))
subplot(3,1,1)
plot(time,aircraft_state_array(:,10),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("Roll Rate (rad/s)");
subplot(3,1,2)
plot(time,aircraft_state_array(:,11),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("Pitch Rate (rad/s)");
subplot(3,1,3)
plot(time,aircraft_state_array(:,12),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("Yaw Rate (rad/s)");

sgtitle("Angular Velocity");

% Control Input Variables
figure(fig(5))
subplot(4,1,1)
plot(time,control_input_array(:,1),col,"Linewidth",1);
grid minor
xlabel("Time (s)");
ylabel("Z Control Input (N)")
subplot(4,1,2)
plot(time,control_input_array(:,2),col,"Linewidth",1);
grid minor
xlabel("Time (s)");
ylabel("L Control Moment (kgm/s^2)")
subplot(4,1,3)
plot(time,control_input_array(:,3),col,"Linewidth",1);
grid minor
xlabel("Time (s)");
ylabel("M Control Moment (kgm/s^2)")
subplot(4,1,4)
plot(time,control_input_array(:,4),col,"Linewidth",1);
grid minor
xlabel("Time (s)");
ylabel("N Control Moment (kgm/s^2)")

% Flight Path
figure(fig(6))
hold on;
plot3(aircraft_state_array(:,1),aircraft_state_array(:,2),aircraft_state_array(:,3),col,'Linewidth',1);
scatter3(aircraft_state_array(1,1),aircraft_state_array(1,2),aircraft_state_array(1,3),'g')
scatter3(aircraft_state_array(end,1),aircraft_state_array(end,2),aircraft_state_array(end,3),'r')
grid minor;
xlabel("X Position (m)");
ylabel('Y Position (m)');
zlabel("Z Position (m)");
title("Aircraft Position");
hold off;

% Figure 7
figure(fig(7))
subplot(8,2,1)
plot(time,aircraft_state_array(:,1),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("x (m)");
subplot(8,2,3)
plot(time,aircraft_state_array(:,2),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("y (m)");
subplot(8,2,5)
plot(time,aircraft_state_array(:,3),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("z (m)");
subplot(8,2,7)
plot(time,aircraft_state_array(:,4),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("\phi (rad)");
subplot(8,2,9)
plot(time,aircraft_state_array(:,5),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("\rho (rad)");
subplot(8,2,11)
plot(time,aircraft_state_array(:,6),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("\theta (rad)");
subplot(8,2,2)
plot(time,aircraft_state_array(:,7),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("uE (m/s)");
subplot(8,2,4)
plot(time,aircraft_state_array(:,8),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("vE (m/s)");
subplot(8,2,6)
plot(time,aircraft_state_array(:,9),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("wE (m/s)");
subplot(8,2,8)
plot(time,aircraft_state_array(:,10),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("p (rad/s)");
subplot(8,2,10)
plot(time,aircraft_state_array(:,11),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("q (rad/s)");
subplot(8,2,12)
plot(time,aircraft_state_array(:,12),col,'Linewidth',1);
grid minor
xlabel('Time (s)');
ylabel("r (rad/s)");
subplot(8,2,13)
plot(time,control_input_array(:,1),col,"linewidth",1);
xlabel("Time");
ylabel("F_1 (N)")
grid minor
subplot(8,2,14)
plot(time,control_input_array(:,3),col,"linewidth",1);
xlabel("Time");
ylabel("F_2 (N)")
grid minor
subplot(8,2,15)
plot(time,control_input_array(:,2),col,"linewidth",1);
xlabel("Time");
ylabel("F_3 (N)")
grid minor
subplot(8,2,16)
plot(time,control_input_array(:,4),col,"linewidth",1);
xlabel("Time");
ylabel("F_4 (N)")
sgtitle("Quadrotor State Variables vs. Time")
grid minor

end