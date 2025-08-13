close all; clear; clc;

const.G = 6.67384*10^-11; %m^3/(kg*s^2)
const.m1 = 5.973*10^24; %kg
const.m2 = 1.899*10^27; %kg

%Works good- according to drawing
x10 = -100000;
y10 = 100000;
v1x0 = 0;
v1y0 = -100000;

x20 = 5000000;
y20 = -1000000;
v2x0 = 10000;
v2y0 = 60000;

% New? Planet flying away
% x10 = -100000;
% y10 = 100000;
% v1x0 = -100000;
% v1y0 = -100000;
% 
% x20 = 5000000;
% y20 = -1000000;
% v2x0 = 10000;
% v2y0 = 60000;


state0 = [x10;y10;x20;y20;v1x0;v1y0;v2x0;v2y0]; 
tspan = [0;200];

f = @(t,state) EOMfunc(t,state,const);

[t,state] = ode45(f,tspan,state0);

lim1 = min(state(:,[1 3]));
lim2 = max(state(:,[1 3]));
lim3 = min(state(:,[2 4]));
lim4 = max(state(:,[2 4]));

lim1 = min(lim1);
lim2 = max(lim2);
lim3 = min(lim3);
lim4 = max(lim4);

lim1 = lim1 + 0.2*lim1;
lim2 = lim2 + 0.1*lim2;
lim3 = lim3 + 0.1*lim3;
lim4 = lim4 + 0.1*lim4;
figure()
hold on;
xlim([lim1 lim2]);
ylim([lim3,lim4]);

for i = 1:size(t)
    scatter(state(i,1),state(i,2),'blue');
    scatter(state(i,3),state(i,4),'red');
    pause(0.001)
end



function [vec] = EOMfunc(~,state,const)
    
    x1 = state(1);
    y1 = state(2);
    x2 = state(3);
    y2 = state(4);
    v1x = state(5);
    v1y = state(6);
    v2x = state(7);
    v2y = state(8);

    r = sqrt((x2-x1)^2+(y2-y1)^2);
    Fg = const.G*const.m1*const.m2/r^2;

    a1x = (1/const.m1)*Fg*(x2-x1)/r;
    a2x = -a1x;

    a1y = (1/const.m1)*Fg*(y2-y1)/r;
    a2y = -a1y;

    vec = [v1x;v1y;v2x;v2y;a1x;a1y;a2x;a2y];


end