function [ddx,xout] = finite_difference(x,y)
    n = length(x);
    ddx = zeros(n,1);
    xout = zeros(n,1);
    for i = 1:n-1
        ddx(i) = (y(i+1)-y(i))/(x(i+1)-x(i));
        xout(i) = (x(i+1)+x(i))/2;
    end
    ddx(n) = ddx(n-1);