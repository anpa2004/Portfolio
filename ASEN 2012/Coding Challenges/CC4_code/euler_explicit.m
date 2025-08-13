function [times,est] = euler_explicit(f,x0,xf,y0,step)
     
     %Calculating N (number of times of iteration)
     N = (xf-x0)/step+1;

     %Creating time vector which outputs from function
     times = linspace(x0,xf,N);
    
    %Setting est(1) equal to the initial y value, preallocating
    est = zeros(length(times),1);
    est(1)=y0;

    %Calculating the estimated value at each moment along step
    for i = 1:(length(times)-1)
        est(i+1) = est(i) + f(times(i),est(i))*step;
    end

end