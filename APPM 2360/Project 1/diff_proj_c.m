clear; close all; clc;

%Time array and constants
tc = 0:.1:24; %s
K = .25; %/s
M_0c1 = 75; %F
M_0c2 = 35; %F

%Integration constant C
Cc = -10+12*K*((cos(-5*pi/12)*(1/K)+(pi/12)*sin(-5*pi/12)*(1/K^2))/(pi^2/(144*K^2)+1));

%Function T(t)
Tc1 = M_0c1 -12*K*((cos((tc-5).*pi./12).*(1/K)+(pi/12)*sin((tc-5).*pi./12).*(1/K^2))./(pi^2/(144*K^2)+1))+Cc./(exp(K*tc));
Mtc = M_0c1 - 12*cos(pi.*(tc-5)./12);
Tc2 = M_0c2 -12*K*((cos((tc-5).*pi./12).*(1/K)+(pi/12)*sin((tc-5).*pi./12).*(1/K^2))./(pi^2/(144*K^2)+1))+Cc./(exp(K*tc));
Mtc2 = M_0c2 - 12*cos(pi.*(tc-5)./12);

%Max value and when
[Tcmax1,tcmax1]=max(Tc1);
[Tcmax2,tcmax2]=max(Tc2);

%Plotting results
figure(1);
subplot(2,1,1);
hold on;
plot(tc,Tc1);
plot(tc,Mtc);
grid on;
xlabel('Time (hrs)');
ylabel('Temperature (^oF)');
title('Temperature change due to varying outside temperature, M_0 = 75^0F');
legend('Temperature of Building (^oF), M_0 = 75^oF','Temperature outside (^oF)');
hold off;

subplot(2,1,2);
hold on;
plot(tc,Tc2);
plot(tc,Mtc2);
grid on;
xlabel('Time (hrs)');
ylabel('Temperature (^oF)');
title('Temperature change due to varying outside temperature, M_0 = 35^0F');
legend('Temperature of Building (^oF), M_0 = 35^oF','Temperature outside (^oF)');
hold off;

%Converting index value to time
tcmaxm = (.1*tcmax1-floor(.1*tcmax1))*60;
tcmaxh = floor(0.1*tcmax1);

tcmaxm2 = (.1*tcmax2-floor(.1*tcmax2))*60;
tcmaxh2 = floor(0.1*tcmax2);

%Printing results
fprintf('(Part C) Max temperature of %2.2f at %g hours and %g minutes for M_0 of 75 F \n',Tcmax1,tcmaxh,tcmaxm);
if Tcmax1<81
    disp('The building does not exceed the max safe temperature of 81 F');
end
if Tcmax1>81
        disp('Equipment will be damaged by the high temperature.');
end

fprintf('(Part C) Max temperature of %2.2f at %g hours and %g minutes for M_0 of 35 F\n',Tcmax2,tcmaxh2,tcmaxm2);
if Tcmax2<81
    disp('The building does not exceed the max safe temperature of 81 F');
end
if Tcmax2>81
        disp('Equipment will be damaged by the high temperature.');
end

