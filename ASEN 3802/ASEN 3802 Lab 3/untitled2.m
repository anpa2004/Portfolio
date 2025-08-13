close all; clear; clc;

x = linspace(0,4,100);
y = sin(x);

[yp,xout] = finite_difference(x,y);

figure()
plot(x,y)
hold on
plot(xout,yp)
