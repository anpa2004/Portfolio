close all; clear; clc;

% x bounds and guesses
a = -1;  
b = 8;
numX = 1000; % Number of x guesses

x_potential = linspace(a,b,numX);

g = @(x) 4*sin(2*x)+3;

Nmax = 1500;
x0 = -1;
tol = 1e-10;

[solution,ier] = fixedpoint(g,x0,tol,Nmax);

function [solution,ier] = fixedpoint(f,x0,tol,Nmax)

    N = 0;
    x1 = f(x0);
    ier = 0;
    while N <= Nmax
        x0 = x1; % Setting x_n to previous x_n+1
        fprintf("n = %g, root = %g\n",N,x1)
        x1 = f(x0);
        N = N+1;

        if abs(x1 - x0) < tol  % Check convergence
            solution = x1;
            return;
        end    
    end

    if abs(x1-x0)> tol
        ier = 1;
        solution = x1;
    end
end