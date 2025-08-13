close all; clear; clc;

%% Constants

%
b = 40;
a0_t = 2*pi;
a0_r = 2*pi;
c_t = 5;
c_r = 5;
aero_t = 0;
aero_r = 0;
geo_t = 5*pi/180;
geo_r = 5*pi/180;

N = 2;

[e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N);

%% PART 2
close all; clear; clc;
a0_t = 2*pi;
a0_r = 2*pi;
c_r = 15;
aero_t = 0;
aero_r = 0;
geo_t = 2*pi/180;
geo_r = 2*pi/180;
N = 10;

ctcr = linspace(0,1,200);
AR = [4,6,8,10];


delta = zeros(length(ctcr),length(AR));
for j = 1:length(AR)
    for i = 1:length(ctcr)
        b = 0.5*AR(j)*(c_r*ctcr(i)+c_r);
        %[e(i,j),~,~,delta(i,j)] = PLLT2(b,a0_t,a0_r,c_r*ctcr(i),c_r,aero_t,aero_r,geo_t,geo_r,N);
        [e(i,j),~,~] = PLLT(b,a0_t,a0_r,c_r*ctcr(i),c_r,aero_t,aero_r,geo_t,geo_r,N);
        delta(i,j) = 1/(e(i,j)) - 1;
    end
end

figure()
hold on
for j = 1:length(AR)
    plot(ctcr,delta(:,j),'Linewidth',1)
end
grid minor
legend("AR = 4","AR = 6","AR = 8","AR = 10",'location','northwest')
xlabel("c_t/c_r")
ylabel("Induced drag factor \delta")
ylim([0,0.18])
title("Aspect Ratio vs. Induced Drag Factor")

figure()
hold on
for j = 1:length(AR)
    plot(ctcr,e(:,j),'Linewidth',1)
end
grid minor
legend("AR = 4","AR = 6","AR = 8","AR = 10",'location','northwest')
xlabel("c_t/c_r")
ylabel("e")
title("Aspect Ratio influence on Spanwise Efficiency ")


