close all; clear; clc;

% Test 
x = linspace(0,2*pi,100);
y = sin(x).*cos(20*x);

X = fft(y);
plot(real(X))
