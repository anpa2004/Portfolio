function [zeta,omegan] = FindDampingNaturalFreq(time,state)

% Taking only the z value of state
x = state(:,7);

% Finding the first peak
[x1,i1] = max(x(150:end));

% Finding the 2nd peak (40 through experimentation)
x2 = max(x(i1+230:end));

% Finding the index when x is x2
i2 = find(x==x2);

% Validation of 2nd peak algorithm
index = linspace(1,length(x),length(x));
figure()
hold on
plot(index(150:end),x(150:end))
plot(index(i2:end),x(i2:end))
yline(x1)
yline(x2)
hold off

%Time difference between the peaks
T = time(i2)-time(i1); 
n = 1; % One period between peaks

% Logarithmic Decriment
delta = 1/n*log(x1/x2);

% Damping Ratio
zeta = delta/sqrt(4*pi^2+delta^2);

% Damped and natural ratios
omegad = 2*pi/T;
omegan = omegad/sqrt(1-zeta^2);
