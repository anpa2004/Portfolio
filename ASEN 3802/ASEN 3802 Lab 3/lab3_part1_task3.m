close all; clear; clc; 


%% Plot cl/alpha for different camber airfoils

n = 20; %Just more than we predicted for the sake of margin
c = 1; %m
alpha = linspace(-5,20,25);

% NACA 0012
digits = [0,0,12];
coords0012 = naca4series(digits,c,n);

% NACA 2412
digits = [2,4,12];
coords2412 = naca4series(digits,c,n);

% NACA 4418
digits = [4,4,12];
coords4412 = naca4series(digits,c,n);

% Finding the cl/alpha graphs
cl0012 = zeros(length(alpha),1);
cl2412 = zeros(length(alpha),1);
cl4412 = zeros(length(alpha),1);
for i = 1:length(alpha)
    cl0012(i) = Vortex_Panel(coords0012.xb,coords0012.yb,alpha(i));
    cl2412(i) = Vortex_Panel(coords2412.xb,coords2412.yb,alpha(i));
    cl4412(i) = Vortex_Panel(coords4412.xb,coords4412.yb,alpha(i));
end

% Plotting Results
figure()
hold on
plot(alpha,cl0012,'Linewidth',1)
plot(alpha,cl2412,'Linewidth',1)
plot(alpha,cl4412,'Linewidth',1)
xline(0)
yline(0)
grid minor
xlabel("\alpha (^\circ)")
ylabel("c_l")
title("c_l vs \alpha")
legend("NACA 0012", "NACA 2412", "NACA 4412",'Location','Northwest')

% Estimating the lift slope and 0 lift AoA

% For lift slop- mostly linear so can use simple slope equation
a0_NACA0012 = (cl0012(3)-cl0012(4))/(alpha(3)-alpha(4));
a0_NACA2412 = (cl2412(3)-cl2412(4))/(alpha(3)-alpha(4));
a0_NACA4412 = (cl4412(3)-cl4412(4))/(alpha(3)-alpha(4));

% Interpolating to find alpha L=0
index_neg = find(cl0012<0,1,"last");
index_pos = find(cl0012>0,1);
alphaL0_NACA0012 = interpinator(cl0012(index_neg),alpha(index_neg),cl0012(index_pos),alpha(index_pos),0);


index_neg = find(cl2412<0,1,"last");
index_pos = find(cl2412>0,1);
alphaL0_NACA2412 = interpinator(cl2412(index_neg),alpha(index_neg),cl2412(index_pos),alpha(index_pos),0);
% Using finite methods
cv = c*ones(length(coords2412.xc),1);
dzdx = finite_difference(coords2412.xc,coords2412.yc);
integrand = 1/pi*dzdx.*(2*coords2412.xc/c).*1/c.*sqrt(2/(coords2412.xc.*(cv-coords2412.xc)));
alphal0_2412 = trapz(integrand(:,10));


index_neg = find(cl4412<0,1,"last");
index_pos = find(cl4412>0,1);
alphaL0_NACA4412 = interpinator(cl4412(index_neg),alpha(index_neg),cl4412(index_pos),alpha(index_pos),0);
cv = c*ones(length(coords4412.xc),1);
dzdx = finite_difference(coords4412.xc,coords4412.yc);
integrand = 1/pi*dzdx.*(2*coords4412.xc/c).*1/c.*sqrt(2/(coords4412.xc.*(cv-coords4412.xc)));
alphal0_4412 = trapz(integrand(:,10));



function [yout] = interpinator(x1,y1,x2,y2,xq)
    m = (y2-y1)/(x2-x1);
    yout = m*(xq-x1)+y1;
end