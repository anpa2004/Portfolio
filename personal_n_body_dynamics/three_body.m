close all; clear; clc;

%Constants
const.m1 = 1.898*10^33; %kg
const.m2 = 5.972*10^24; %kg
const.m3 = 3.4810^22; %kg
const.G = 6.67384*10^-11; %m^3/(kg*s^2)

x10 = -100000;
y10 = 100000;
z10 = 10020;
v1x0 = 0;
v1y0 = -100000;
v1z0 = -10000;

x20 = 5000000;
y20 = -1000000;
z20 = -100403;
v2x0 = 10000;
v2y0 = 60000;
v2z0 = 10000;

x30 = -1002000;
y30 = 50002005;
z30 = -423401;
v3x0 = 60000;
v3y0 = -102001;
v3z0 = 0;

state0 = [x10;y10;z0;x20;y20;z20;x30;y30;z30;v1x0;v1y0;v1z0;v2x0;v2y0;v2z0;v3x0;v3y0;v3z0];






function [vec] = EOMfunc(t,state,const)
    x1 = state(1);
    y1 = state(2);
    z1 = state(3);
    x2 = state(4);
    y2 = state(5);
    z2 = state(6);
    x3 = state(7);
    y3 = state(8);
    z3 = state(9);
    v1x = state(10);
    v1y = state(11);
    v1z = state(12);
    v2x = state(13);
    v2y = state(14);
    v2z = state(15);
    v3x = state(16);
    v3y = state(17);
    v3z = state(18);


    r12 = sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);
    r13 = sqrt((x3-x1)^2+(y3-y1)^2+(z3-z1)^2);
    r23 = sqrt((x2-x3)^2+(y2-y3)^2+(z3-z2)^2);

    Fg12 = const.m1*const.m2*const.G/(r12^2);
    Fg13 = const.m1*const.m3*const.G/(r13^2);
    Fg23 = const.m2*const.m3*const.G/(r23^2);

    Fg12x = Fg12*(x2-x1)/r12;
    Fg12y = Fg12*(y2-y1)/r12;
    Fg12z = Fg12*(z2-z1)/r12;

    Fg13x = Fg13*(x3-x1)/r13;
    Fg13y = Fg13*(y3-y1)/r13;
    Fg13z = Fg13*(z3-z1)/r13;

    Fg23x = Fg23*(x3-x2)/r23;
    Fg23y = Fg23*(y3-y2)/r23;
    Fg23z = Fg23*(z3-z2)/r23;

    a1x = (1/const.m1)*Fg*(x2-x1)/r;
    a2x = -a1x;



end