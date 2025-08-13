close all; clear; clc;

%% Givens 
b = 35 + 10/12; %ft
c_r = 5 + 4/12; %ft
c_t = 7 + 7/12; %ft
a0_t = 0.1185 * 180/pi;
a0_r = 0.1186 * 180/pi; % FROM PART 1 OF LAB
aero_t = 0;
aero_r = 0;
nums_root = [2,4,12];
nums_tip = [0,0,12];
geo_root = 0;
geo_tip = geo_root + 2*pi/180; %rad
N = 20;

%% PART 1: Coefficient of lift
airfoil_root = naca4series(nums_root,c_r,20);
airfoil_tip = naca4series(nums_tip,c_t,20);

alpha = linspace(-5,10,100);

for i = 1:length(alpha)
    cl_tip(i) = Vortex_Panel(airfoil_tip.xb,airfoil_tip.yb,alpha(i)+2);
    cl_root(i) = Vortex_Panel(airfoil_root.xb,airfoil_root.yb,alpha(i));
    [~,c_L(i),c_Di(i)] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,(alpha(i) + 2)*pi/180,alpha(i)*pi/180,N);
end

%Finding a0
lift_slope = (c_L(3)-c_L(2))/(alpha(3)-alpha(2))

% Plotting results
figure()
hold on;
grid minor
plot(alpha,c_L,'Linewidth',1)
xline(0)
yline(0)
title("c_l vs \alpha for Whole Aircraft")

figure()
hold on;
grid minor 
plot(alpha,c_Di,"Linewidth",1)
xlabel("\alpha (deg)")
ylabel("c_d")
title("Induced Drag versus Alpha (From PLLT)")
