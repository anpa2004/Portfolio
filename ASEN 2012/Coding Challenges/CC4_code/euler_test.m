close all; clear; clc;

x0=0;
y0=0;
xf = 0.5;

dydt = @(t,y) cos(t);

[tout,yout] = euler_explicit(0.1,x0,xf,y0,dydt);

