%% ASEN 3802 Lab 2 Task 1

clc;clear;close all;

const.T0 = 12; % [C]
const.H_an = 98.6;  % [C/m]
const.k = 130; % [W/mK]
const.rho = 2810; % [kg/m^3]
const.c_p = 960; % [J/kgK]
const.alpha = const.k/(const.rho*const.c_p); % [m^2/s]
const.x = 0.1238; % [m]
const.L = const.x+0.0508; % [m]

n = 10;
t1 = 1; % [s]
t2 = 1000; % [s]

u_1s = heatdistribution(t1,n,const);
u_1000s = heatdistribution(t2,n,const);

Fo_1s = (const.alpha*t1)/const.L^2;
Fo_1000s = (const.alpha*t2)/const.L^2;

figure 
hold on
plot(0:n,u_1s,'linewidth',1)
plot(0:n,u_1000s,'linewidth',1)
grid minor
title('Temperature at time t vs n')
xlabel('n')
ylabel('Temperature at t [C]')
legend('t = 1s','t = 1000s')


function u = heatdistribution(t,n,const)
    sum = 0;
    u(1) = const.T0+const.H_an*const.x;
    for i = 1:n
        b_n = (8*const.H_an*const.L*(-1)^i)/((2*i-1)^2*pi^2);
        lambda_n = pi*(2*i-1)/(2*const.L);
        sum = sum + b_n*sin(lambda_n*const.x)*exp(-1*lambda_n^2*const.alpha*t);
        u(i+1) = const.T0+const.H_an*const.x + sum;
    end
end
