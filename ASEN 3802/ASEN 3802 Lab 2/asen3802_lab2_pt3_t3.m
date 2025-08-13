clc;clear;close all;

al26 = readmatrix("Aluminum_26V_250mA");
al28 = readmatrix("Aluminum_28V_269mA");
br26 = readmatrix("Brass_26V_245mA");
br29 = readmatrix("Brass_29V_273mA");
st = readmatrix("Steel_21V_192mA");

%% Material and Geometric Properties
dx = 0.0127; %m
x0 = 0.034925;
xpos = x0:dx:x0 + 0.1016;
% breaking data into channels
% cell row is which thermocouple
al26_c = cell(8,1);
al28_c = cell(8,1);
br26_c = cell(8,1);
br29_c = cell(8,1);
st_c = cell(8,1);
al26_t = al26(:,1);
al28_t = al28(:,1);
br26_t = br26(:,1);
br29_t = br29(:,1);
st_t = st(:,1);
for i = 1:8
al26_c{i} = al26(:,i+1);
al28_c{i} = al28(:,i+1);
br26_c{i} = br26(:,i+1);
br29_c{i} = br29(:,i+1);
st_c{i} = st(:,i+1);
end

al26_notime = al26(:,2:end);
al28_notime = al28(:,2:end);
br26_notime = br26(:,2:end);
br29_notime = br29(:,2:end);
st_notime = st(:,2:end);

x1 = 0.0349; % [m]
x2 = 0.0476; % [m]
x3 = 0.0603; % [m]
x4 = 0.073; % [m]
x5 = 0.0857; % [m]
x6 = 0.0984; % [m]
x7 = 0.1111; % [m]
x8 = 0.1238; % [m]
x = [x1,x2,x3,x4,x5,x6,x7,x8];

const_st.H_exp = 277; % [C/m]
const_st.H_an = 491.2; % [C/m]
const_st.M_exp = 35.2; % [C/m]
const_st.T0 = 11; % [C]
const_st.k = 16.2; % [W/mK]
const_st.cp = 500; % [J/kgK]
const_st.rho = 8000; % [kg/m^3]
const_st.alpha = const_st.k/(const_st.rho*const_st.cp); % [m^2/s]
const_st.L = x8+0.0508; % [m]
const_st.alpha2 = 7.36e-6;

const_br26.H_exp = 114; % [C/m]
const_br26.H_an = 109.3; % [C/m]
const_br26.M_exp = 7; % [C/m]
const_br26.T0 = 12; % [C]
const_br26.k = 115; % [W/mK]
const_br26.cp = 380; % [J/kgK]
const_br26.rho = 8500; % [kg/m^3]
const_br26.alpha = const_br26.k/(const_br26.rho*const_br26.cp); % [m^2/s]
const_br26.L = x8+0.0508; % [m]
const_br26.alpha2 = 2.3e-5;

 
const_br29.H_exp = 139.8; % [C/m]
const_br29.H_an = 135.9; % [C/m]
const_br29.M_exp = 6.6; % [C/m]
const_br29.T0 = 12; % [C]
const_br29.k = 115; % [W/mK]
const_br29.cp = 380; % [J/kgK]
const_br29.rho = 8500; % [kg/m^3]
const_br29.alpha = const_br29.k/(const_br29.rho*const_br29.cp); % [m^2/s]
const_br29.L = x8+0.0508; % [m]
const_br29.alpha2 = 2.3e-5;

const_al26.H_exp = 54.3; % [C/m]
const_al26.H_an = 98.6; % [C/m]
const_al26.M_exp = -3.4; % [C/m]
const_al26.T0 = 12; % [C]
const_al26.k = 130; % [W/mK]
const_al26.cp = 960; % [J/kgK]
const_al26.rho = 2810; % [kg/m^3]
const_al26.alpha = const_al26.k/(const_al26.rho*const_al26.cp); % [m^2/s]
const_al26.L = x8+0.0508; % [m]
const_al26.alpha2 = 3.6e-5;

const_al28.H_exp = 66.8; % [C/m]
const_al28.H_an = 114.3; % [C/m]
const_al28.M_exp = 0.3; % [C/m]
const_al28.T0 = 12; % [C]
const_al28.k = 130; % [W/mK]
const_al28.cp = 960; % [J/kgK]
const_al28.rho = 2810; % [kg/m^3]
const_al28.alpha = const_al28.k/(const_al28.rho*const_al28.cp); % [m^2/s]
const_al28.L = x8+0.0508; % [m]
const_al28.alpha2 = 4.97e-5;

alpha_br26 = linspace(0,2*const_br26.alpha,100);
alpha_br29 = linspace(0,2*const_br29.alpha,100);
alpha_al26 = linspace(0,2*const_al26.alpha,100);
alpha_al28 = linspace(0,2*const_al28.alpha,100);
alpha_st = linspace(0,2*const_st.alpha,100);

t_st = 0:10:12510; % [s]
t_al26 = 0:10:1670; % [s]
t_al28 = 0:10:2670; % [s]
t_br26 = 0:10:4940; % [s]
t_br29 = 0:10:4930; % [s]

%% Finding Analytical Solutions

% Analytical solutions for model IA
for i = 1:length(t_st)
    u_stIA(i,:) = heatdistributionIA(t_st(i),x8,const_st);
end
for i = 1:length(t_al26)
    u_al26IA(i,:) = heatdistributionIA(t_al26(i),x8,const_al26);
end
for i = 1:length(t_al28)
    u_al28IA(i,:) = heatdistributionIA(t_al28(i),x8,const_al28);
end
for i = 1:length(t_br26)
    u_br26IA(i,:) = heatdistributionIA(t_br26(i),x8,const_br26);
end
for i = 1:length(t_br29)
    u_br29IA(i,:) = heatdistributionIA(t_br29(i),x8,const_br29);
end

% Analytical solutions for model IB
for i = 1:length(t_st)
    u_stIB(i,:) = heatdistributionIB(t_st(i),x8,const_st);
end
for i = 1:length(t_al26)
    u_al26IB(i,:) = heatdistributionIB(t_al26(i),x8,const_al26);
end
for i = 1:length(t_al28)
    u_al28IB(i,:) = heatdistributionIB(t_al28(i),x8,const_al28);
end
for i = 1:length(t_br26)
    u_br26IB(i,:) = heatdistributionIB(t_br26(i),x8,const_br26);
end
for i = 1:length(t_br29)
    u_br29IB(i,:) = heatdistributionIB(t_br29(i),x8,const_br29);
end

% Analytical Solutions for model II
for i = 1:length(t_st)
    u_stII(i,:) = heatdistributionII(t_st(i),x8,const_st);
end
for i = 1:length(t_al26)
    u_al26II(i,:) = heatdistributionII(t_al26(i),x8,const_al26);
end
for i = 1:length(t_al28)
    u_al28II(i,:) = heatdistributionII(t_al28(i),x8,const_al28);
end
for i = 1:length(t_br26)
    u_br26II(i,:) = heatdistributionII(t_br26(i),x8,const_br26);
end
for i = 1:length(t_br29)
    u_br29II(i,:) = heatdistributionII(t_br29(i),x8,const_br29);
end

% Analytical Solution for model III
for i = 1:length(t_st)
    u_stIII(i,:) = heatdistributionIII(t_st(i),x8,const_st);
end
for i = 1:length(t_al26)
    u_al26III(i,:) = heatdistributionIII(t_al26(i),x8,const_al26);
end
for i = 1:length(t_al28)
    u_al28III(i,:) = heatdistributionIII(t_al28(i),x8,const_al28);
end
for i = 1:length(t_br26)
    u_br26III(i,:) = heatdistributionIII(t_br26(i),x8,const_br26);
end
for i = 1:length(t_br29)
    u_br29III(i,:) = heatdistributionIII(t_br29(i),x8,const_br29);
end


%% Data analysis
dx = 0.0127; %m
x0 = 0.034925;
xpos = x0:dx:x0 + 0.1016;
% breaking data into channels
% cell row is which thermocouple
al26_c = cell(8,1);
al28_c = cell(8,1);
br26_c = cell(8,1);
br29_c = cell(8,1);
st_c = cell(8,1);
al26_t = al26(:,1);
al28_t = al28(:,1);
br26_t = br26(:,1);
br29_t = br29(:,1);
st_t = st(:,1);
for i = 1:8
al26_c{i} = al26(:,i+1);
al28_c{i} = al28(:,i+1);
br26_c{i} = br26(:,i+1);
br29_c{i} = br29(:,i+1);
st_c{i} = st(:,i+1);
end

%% Plotting Results\

% AL 26
plotvec = al26_t(1:25:end);
err = 2*ones(length(plotvec),1);

figure()
hold on
p1=plot(al26_t,al26_c{8},'LineWidth',1);
p2=errorbar(al26_t(1:25:end),al26_c{8}(1:25:end),err,'k');
p3=plot(t_al26,u_al26IA(:,8),'LineWidth',1);
p4=plot(t_al26,u_al26IB(:,8),'LineWidth',1);
p5=plot(t_al26,u_al26II(:,8),'LineWidth',1);
p6=plot(t_al26,u_al26III(:,8),'LineWidth',1);
grid minor
xlabel("Time (s)")
ylabel('Temperature (^\circC)')
legend([p1,p2,p3,p4,p5,p6],{'Experimental Results','\pm2 ^\circC error','Model IA','Model 1B','Model II','Model III'},'Location','southeast')
legend
title("Model Comparison: Aluminum - 26V,250mA")


% Steel
plotvec = st_t(1:25:end);
err = 2*ones(length(plotvec),1);

figure()
hold on
p1=plot(st_t,st_c{8},'LineWidth',1);
p2=errorbar(st_t(1:25:end),st_c{8}(1:25:end),err,'k');
p3=plot(t_st,u_stIA(:,8),'LineWidth',1);
p4=plot(t_st,u_stIB(:,8),'LineWidth',1);
p5=plot(t_st,u_stII(:,8),'LineWidth',1);
p6=plot(t_st,u_stIII(:,8),'LineWidth',1);
grid minor
xlabel("Time (s)")
ylabel('Temperature (^\circC)')
legend([p1,p2,p3,p4,p5,p6],{'Experimental Results','\pm2 ^\circC error','Model IA','Model 1B','Model II','Model III'},'Location','southeast')
legend
title("Model Comparison: Steel - 21V,192mA")


% AL 28
plotvec = al28_t(1:25:end);
err = 2*ones(length(plotvec),1);

figure()
hold on
p1=plot(al28_t,al28_c{8},'LineWidth',1);
p2=errorbar(al28_t(1:25:end),al28_c{8}(1:25:end),err,'k');
p3=plot(t_al28,u_al28IA(:,8),'LineWidth',1);
p4=plot(t_al28,u_al28IB(:,8),'LineWidth',1);
p5=plot(t_al28,u_al28II(:,8),'LineWidth',1);
p6=plot(t_al28,u_al28III(:,8),'LineWidth',1);
grid minor
xlabel("Time (s)")
ylabel('Temperature (^\circC)')
legend([p1,p2,p3,p4,p5,p6],{'Experimental Results','\pm2 ^\circC error','Model IA','Model 1B','Model II','Model III'},'Location','southeast')
legend
title("Model Comparison: Aluminum - 28V,269mA")

% Br 26
plotvec = br26_t(1:25:end);
err = 2*ones(length(plotvec),1);

figure()
hold on
p1=plot(br26_t,br26_c{8},'LineWidth',1);
p2=errorbar(br26_t(1:25:end),br26_c{8}(1:25:end),err,'k');
p3=plot(t_br26,u_br26IA(:,8),'LineWidth',1);
p4=plot(t_br26,u_br26IB(:,8),'LineWidth',1);
p5=plot(t_br26,u_br26II(:,8),'LineWidth',1);
p6=plot(t_br26,u_br26III(:,8),'LineWidth',1);
grid minor
xlabel("Time (s)")
ylabel('Temperature (^\circC)')
legend([p1,p2,p3,p4,p5,p6],{'Experimental Results','\pm2 ^\circC error','Model IA','Model 1B','Model II','Model III'},'Location','southeast')
legend
title("Model Comparison: Brass - 26V,245mA")

% Br 28
plotvec = br29_t(1:25:end);
err = 2*ones(length(plotvec),1);

figure()
hold on
p1=plot(br29_t,br29_c{8},'LineWidth',1);
p2=errorbar(br29_t(1:25:end),br29_c{8}(1:25:end),err,'k');
p3=plot(t_br29,u_br29IA(:,8),'LineWidth',1);
p4=plot(t_br29,u_br29IB(:,8),'LineWidth',1);
p5=plot(t_br29,u_br29II(:,8),'LineWidth',1);
p6=plot(t_br29,u_br29III(:,8),'LineWidth',1);
grid minor
xlabel("Time (s)")
ylabel('Temperature (^\circC)')
legend([p1,p2,p3,p4,p5,p6],{'Experimental Results','\pm2 ^\circC error','Model IA','Model 1B','Model II','Model III'},'Location','southeast')
legend
title("Model Comparison: Brass - 29V,273mA")

%% Functions
function u = heatdistributionIA(t,x,const)
    sum = [0,0,0,0,0,0,0,0];
    for n = 1:10
        b_n = (8*const.H_an*const.L*(-1)^n)/((2*n-1)^2*pi^2);
        lambda_n = pi*(2*n-1)/(2*const.L);
        sum = sum + b_n.*sin(lambda_n.*x).*exp(-1*lambda_n^2*const.alpha*t);
    end
    u = const.T0+const.H_an.*x + sum;
end

function u = heatdistributionIB(t,x,const)
    sum = [0,0,0,0,0,0,0,0];
    for n = 1:10
        b_n = (8*const.H_exp*const.L*(-1)^n)/((2*n-1)^2*pi^2);
        lambda_n = pi*(2*n-1)/(2*const.L);
        sum = sum + b_n.*sin(lambda_n.*x).*exp(-1*lambda_n^2*const.alpha*t);
    end
    u = const.T0+const.H_exp.*x + sum;
end

function u = heatdistributionII(t,x,const)
    sum = [0,0,0,0,0,0,0,0];
    for n = 1:10
        b_n = (8*(const.M_exp-const.H_exp)*const.L*(-1)^(n+1))/((2*n-1)^2*pi^2);
        lambda_n = pi*(2*n-1)/(2*const.L);
        sum = sum + b_n.*sin(lambda_n.*x).*exp(-1*lambda_n^2*const.alpha*t);
    end
    u = const.T0+const.H_exp.*x + sum;
end

function u = heatdistributionIII(t,x,const)
     sum = [0,0,0,0,0,0,0,0];
    for n = 1:10
        b_n = (8*const.H_exp*const.L*(-1)^n)/((2*n-1)^2*pi^2);
        lambda_n = pi*(2*n-1)/(2*const.L);
        sum = sum + b_n.*sin(lambda_n.*x).*exp(-1*lambda_n^2*const.alpha2*t);
    end
    u = const.T0+const.H_exp.*x + sum;
end
