close all; clear; clc;

c = 2;%m
n = 50; %num panels + 1

%% NACA 0021
digits = [0,0,21];

coords = naca4series(digits,c,n);
coords0021.x = cat(1,flip(coords.xL),coords.xU);
coords0021.y = cat(1,flip(coords.yL),coords.yU);

figure()
hold on
plot(coords0021.x,coords0021.y,'Linewidth',1)
plot(coords.xc,coords.yc,'Linewidth',1)
legend("Airfoil surface","Camber line")
grid minor
xlabel("Length along chord (m)")
ylabel("Height (m)")
title("NACA 0021 Airfoil")
ymin = min(coords.xL)-0.4;
ylim([ymin,ymin+c ])
xlim([0,c])

%% NACA 2421
digits = [2,4,21];

coords = naca4series(digits,c,n);
coords2421.x = cat(1,flip(coords.xL),coords.xU);
coords2421.y = cat(1,flip(coords.yL),coords.yU);

figure()
hold on
plot(coords2421.x,coords2421.y,'Linewidth',1)
plot(coords.xc,coords.yc,'Linewidth',1)
legend("Airfoil surface","Camber line")
grid minor
xlabel("Length along chord (m)")
ylabel("Height (m)")
title("NACA 2421 Airfoil")
ymin = min(coords.xL)-0.4;
ylim([ymin,ymin+c ])
xlim([0,c])




%%
digits = [8,2,0];

coords = naca4series(digits,c,n);
coords2421.x = cat(1,flip(coords.xL),coords.xU);
coords2421.y = cat(1,flip(coords.yL),coords.yU);

figure()
hold on
plot(coords2421.x,coords2421.y,'Linewidth',1)
plot(coords.xc,coords.yc,'Linewidth',1)
legend("Airfoil surface","Camber line")
grid minor
xlabel("Length along chord (m)")
ylabel("Height (m)")
title("NACA 2421 Airfoil")
ymin = min(coords.xL)-0.4;
axis equal




