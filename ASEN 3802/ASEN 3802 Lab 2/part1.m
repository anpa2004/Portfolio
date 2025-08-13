close all; clear; clc;

%Storing Data
al1 = readmatrix("Aluminum_26V_250mA.txt");
al2 = readmatrix("Aluminum_28V_269mA.txt");
br1 = readmatrix("Brass_26V_245mA.txt");
br2 = readmatrix("Brass_29V_273mA.txt");
st1 = readmatrix("Steel_21V_192mA.txt");

%% Properties:
% Thermocouple locations
dx = 0.0127; %m
x0 = 0.034925; 
xpos = x0:dx:x0 + 0.1016 - dx;

% Conductivity
k_al = 130;
k_br = 115;
k_st = 16.2;

% Cross secitonal area
A = pi*(0.0127)^2;

% Heat transfer
Q_al1 = 26*0.25;
Q_al2 = 28*0.269;
Q_br1 = 26*0.245;
Q_br2 = 29*0.273;
Q_st1 = 21*0.192; 

% Steady State slope (analytical)
H_al1 = Q_al1/(k_al*A);
H_al2 = Q_al2/(k_al*A);
H_br1 = Q_br1/(k_br*A);
H_br2 = Q_br2/(k_br*A);
H_st1 = Q_st1/(k_st*A);

%% Experimental Data
% breaking data into channels
% cell row is which thermocouple
%prealocation
al1_c = cell(8,1);
al2_c = cell(8,1);
br1_c = cell(8,1);
br2_c = cell(8,1);
st1_c = cell(8,1);
al1_ss_start = zeros(8,1);
al2_ss_start = zeros(8,1);
br1_ss_start = zeros(8,1);
br2_ss_start = zeros(8,1);
st1_ss_start = zeros(8,1);
al1_ss_end = zeros(8,1);
al2_ss_end = zeros(8,1);
br1_ss_end = zeros(8,1);
br2_ss_end = zeros(8,1);
st1_ss_end = zeros(8,1);
al1_t = al1(:,1);
al2_t = al2(:,1);
br1_t = br1(:,1);
br2_t = br2(:,1);
st1_t = st1(:,1);

for i = 1:8
    % Storing experimental data
    al1_c{i} = al1(:,i+1);
    al2_c{i} = al2(:,i+1);
    br1_c{i} = br1(:,i+1);
    br2_c{i} = br2(:,i+1);
    st1_c{i} = st1(:,i+1);

    % Creating initial steady state vector
    al1_ss_start(i) = al1_c{i}(1);
    al2_ss_start(i) = al2_c{i}(1);
    br1_ss_start(i) = br1_c{i}(1);
    br2_ss_start(i) = br2_c{i}(1);
    st1_ss_start(i) = st1_c{i}(1);

    % Creating Final Steady State vector
    al1_ss_end(i) = al1_c{i}(end);
    al2_ss_end(i) = al2_c{i}(end);
    br1_ss_end(i) = br1_c{i}(end);
    br2_ss_end(i) = br2_c{i}(end);
    st1_ss_end(i) = st1_c{i}(end);
end

%% Plotting the Data
figure()
hold on
for i = 1:8
    plot(al1_t,al1_c{i})
end
grid minor
xlabel("Time (s)")
ylabel('Temperature (^\circ C)')
legend('Channel 1','Channel 2','Channel 3','Channel 4','Channel 5','Channel 6','Channel 7','Channel 8','Location','northwest')
title("Aluminum - 26V,250mA")

figure()
hold on
for i = 1:8
    plot(al2_t,al2_c{i})
end
grid minor
xlabel("Time (s)")
ylabel('Temperature (^\circ C)')
legend('Channel 1','Channel 2','Channel 3','Channel 4','Channel 5','Channel 6','Channel 7','Channel 8','Location','northwest')
title("Aluminum - 28V,269mA")

figure()
hold on
for i = 1:8
    plot(br1_t,br1_c{i})
end
grid minor
xlabel("Time (s)")
ylabel('Temperature (^\circ C)')
legend('Channel 1','Channel 2','Channel 3','Channel 4','Channel 5','Channel 6','Channel 7','Channel 8','Location','northwest')
title("Brass - 26V,245mA")

figure()
hold on
for i = 1:8
    plot(br2_t,br2_c{i})
end
grid minor
xlabel("Time (s)")
ylabel('Temperature (^\circ C)')
legend('Channel 1','Channel 2','Channel 3','Channel 4','Channel 5','Channel 6','Channel 7','Channel 8','Location','northwest')
title("Aluminum - 29V,273mA")

figure()
hold on
for i = 1:8
    plot(st1_t,st1_c{i})
end
grid minor
xlabel("Time (s)")
ylabel('Temperature (^\circ C)')
legend('Channel 1','Channel 2','Channel 3','Channel 4','Channel 5','Channel 6','Channel 7','Channel 8','Location','northwest')
title("Steel - 21V,192mA")


%% Creating Linear Fits of Steady States

% Generating Fits
x_query = linspace(0,0.0254+xpos(end),10);

% Start steady state
p1s = polyfit(xpos,al1_ss_start,1);
al1_fits = polyval(p1s,x_query);
p2s = polyfit(xpos,al2_ss_start,1);
al2_fits = polyval(p2s,x_query);
p3s = polyfit(xpos,br1_ss_start,1);
br1_fits = polyval(p3s,x_query);
p4s = polyfit(xpos,br2_ss_start,1);
br2_fits = polyval(p4s,x_query);
p5s = polyfit(xpos,st1_ss_start,1);
st1_fits = polyval(p5s,x_query);

% End Steady State
% 1st Aluminum sample
p1e = polyfit(xpos,al1_ss_end,1);
al1_fite = polyval(p1e,x_query);
al1_T0 = p1e(2);

% Second AL sample
p2e = polyfit(xpos,al2_ss_end,1);
al2_fite = polyval(p2e,x_query);
al2_T0 = p2e(2);

% 1st Brass sample
p3e = polyfit(xpos,br1_ss_end,1);
br1_fite = polyval(p3e,x_query);
br1_T0 = p3e(2);

% 2nd Brass Sample
p4e = polyfit(xpos,br2_ss_end,1);
br2_fite = polyval(p4e,x_query);
br2_T0 = p4e(2);

% Steel Sample
p5e = polyfit(xpos,st1_ss_end,1);
st1_fite = polyval(p5e,x_query);
st1_T0 = p5e(2);


%% Plots of start steady state
figure()
hold on
scatter(xpos,al1_ss_start)
plot(x_query,al1_fits,'LineWidth',1)
legend("Data","Fit")
grid minor
xlabel("x Position (m)")
ylabel("Temperature (^\circC)")
title("Initial Steady State: Aluminum - 26V,250mA")

figure()
hold on
scatter(xpos,al2_ss_start)
plot(x_query,al2_fits,'LineWidth',1)
legend("Data","Fit")
grid minor
xlabel("x Position (m)")
ylabel("Temperature (^\circC)")
title("Initial Steady State: Aluminum - 28V,269mA")

figure()
hold on
scatter(xpos,br1_ss_start)
plot(x_query,br1_fits,'LineWidth',1)
legend("Data","Fit")
grid minor
xlabel("x Position (m)")
ylabel("Temperature (^\circC)")
title("Initial Steady State: Brass- 26V,245mA")

figure()
hold on
scatter(xpos,br2_ss_start)
plot(x_query,br2_fits,'LineWidth',1)
legend("Data","Fit")
grid minor
xlabel("x Position (m)")
ylabel("Temperature (^\circC)")
title("Initial Steady State: Brass- 29V,273mA")


figure()
hold on
scatter(xpos,st1_ss_start)
plot(x_query,st1_fits,'LineWidth',1)
legend("Data","Fit")
grid minor
xlabel("x Position (m)")
ylabel("Temperature (^\circC)")
title("Initial Steady State: Steel- 21V,192mA")


%% Plots of end steady state

% Calculating analytical solution
al1e_an = al1_T0 + H_al1.*x_query;
al2e_an = al2_T0 + H_al2.*x_query;
br1e_an = br1_T0 + H_br1.*x_query;
br2e_an = br2_T0 + H_br2.*x_query;
st1e_an = st1_T0 + H_st1.*x_query;


figure()
hold on
scatter(xpos,al1_ss_end)
plot(x_query,al1_fite,'LineWidth',1)
plot(x_query,al1e_an,'LineWidth',1)
legend("Data","Fit","Analytical")
grid minor
xlabel("x Position (m)")
ylabel("Temperature (^\circC)")
title("Final Steady State: Aluminum - 26V,250mA")

figure()
hold on
scatter(xpos,al2_ss_end)
plot(x_query,al2_fite,'LineWidth',1)
plot(x_query,al2e_an,'LineWidth',1)
legend("Data","Fit","Analytical")
grid minor
xlabel("x Position (m)")
ylabel("Temperature (^\circC)")
title("Final Steady State: Aluminum - 28V,269mA")

figure()
hold on
scatter(xpos,br1_ss_end)
plot(x_query,br1_fite,'LineWidth',1)
plot(x_query,br1e_an,'LineWidth',1)
legend("Data","Fit","Analytical")
grid minor
xlabel("x Position (m)")
ylabel("Temperature (^\circC)")
title("Final Steady State: Brass- 26V,245mA")

figure()
hold on
scatter(xpos,br2_ss_end)
plot(x_query,br2_fite,'LineWidth',1)
plot(x_query,br2e_an,'LineWidth',1)
legend("Data","Fit","Analytical")
grid minor
xlabel("x Position (m)")
ylabel("Temperature (^\circC)")
title("Final Steady State: Brass- 29V,273mA")


figure()
hold on
scatter(xpos,st1_ss_end)
plot(x_query,st1_fite,'LineWidth',1)
plot(x_query,st1e_an,'LineWidth',1)
legend("Data","Fit","Analytical")
grid minor
xlabel("x Position (m)")
ylabel("Temperature (^\circC)")
title("Final Steady State: Steel- 21V,192mA")



