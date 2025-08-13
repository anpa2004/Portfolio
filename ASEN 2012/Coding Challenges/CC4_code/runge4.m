function[x,y] = runge4(f,x0,xf,y0,delta)
    
    %N is number of points or how many times rk4 iterates
    N = (xf-x0)/delta;
    
    %Creating the output vector for the x axis
    x = x0:delta:xf;

    %Preallocating for efficiency
    y = zeros(N,1);

    %The first entry in y matrix is y0
    y(1) = y0;

    %iterating every time through the steps
    for i = 1:(N)
        y1 = f(x(i),y(i));
        y2 = f((x(i) + 0.5 * delta),(y(i) + 0.5 * y1 * delta));
        y3 = f((x(i) + 0.5 * delta),(y(i) + 0.5 * y2 * delta));
        y4 = f((x(i) + delta),(y(i) + delta * y3));

        y(i+1) = y(i) + (1/6)*(y1+ 2*y2 + 2*y3 + y4)*delta;
    end


end