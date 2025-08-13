close all; clear; clc; 

%% TASK 2: Convergence Study

% Properties
n = linspace(5,500,200);
c = 1;
digits = [0,0,6];
alpha = 10; %deg

% Calculating the cl for every n value
for i = 1:length(n)
    coords0006 = naca4series(digits,c,n(i));
    cl0006(i) = Vortex_Panel(coords0006.xb,coords0006.yb,alpha);
end

% Calculating the relative error
error = zeros(length(n)-1,1);
for i = 1:length(n)-1
    error(i) = (cl0006(i)-cl0006(end))/cl0006(end);
end

% Plotting Results
figure()
hold on
plot(n,cl0006,"LineWidth",1)
yline(0.99*cl0006(end))
yline(1.01*cl0006(end))
legend("Data","\pm 1% ","")
grid minor
xlabel("Number of Boundary points")
ylabel("C_l")
title("Coefficient of Lift vs Number of Panels for NACA 0006")


index = find(cl0006>cl0006(end)*.99,1);

%% Plot cl/alpha for different thickness airfoils
clear; close all

n = 20; %Just more than we predicted for the sake of margin
c = 1; %m
alpha = linspace(-5,20,25);

% NACA 0006
digits = [0,0,6];
coords0006 = naca4series(digits,c,n);

% NACA 0012
digits = [0,0,12];
coords0012 = naca4series(digits,c,n);

% NACA 0018
digits = [0,0,18];
coords0018 = naca4series(digits,c,n);

% Finding the cl/alpha graphs
cl0006 = zeros(length(alpha),1);
cl0012 = zeros(length(alpha),1);
cl0018 = zeros(length(alpha),1);
for i = 1:length(alpha)
    cl0006(i) = Vortex_Panel(coords0006.xb,coords0006.yb,alpha(i));
    cl0012(i) = Vortex_Panel(coords0012.xb,coords0012.yb,alpha(i));
    cl0018(i) = Vortex_Panel(coords0018.xb,coords0018.yb,alpha(i));
end

% Plotting Results
figure()
hold on
plot(alpha,cl0006,'Linewidth',1)
plot(alpha,cl0012,'Linewidth',1)
plot(alpha,cl0018,'Linewidth',1)
xline(0)
yline(0)
grid minor
xlabel("\alpha (^\circ)")
ylabel("c_l")
title("c_l vs \alpha")
legend("NACA 0006", "NACA 0012", "NACA 0018",'Location','Northwest')

% Estimating the lift slope and 0 lift AoA

% For lift slop- mostly linear so can use simple slope equation
a0_NACA0006 = (cl0006(3)-cl0006(4))/(alpha(3)-alpha(4));
a0_NACA0012 = (cl0012(3)-cl0012(4))/(alpha(3)-alpha(4));
a0_NACA0018 = (cl0018(3)-cl0018(4))/(alpha(3)-alpha(4));
    

%% Estimating alphaL=0
cv = c*ones(length(coords0006.xc),1);
dzdx = finite_difference(coords0006.xc,coords0006.yc);
integrand = -1/pi*dzdx.*(2*coords0006.xc/c).*1/c.*sqrt(2/(coords0006.xc.*(cv-coords0006.xc)));
alphal0_0006 = trapz(integrand(:,1));

