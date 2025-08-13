%% ASEN 3802 Lab 2 Part 3 Task 1

clc;clear;close all;

al26 = readmatrix("Aluminum_26V_250mA");
al28 = readmatrix("Aluminum_28V_269mA");
br26 = readmatrix("Brass_26V_245mA");
br29 = readmatrix("Brass_29V_273mA");
st = readmatrix("Steel_21V_192mA");

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
const_st.T0 = 11; % [C]
const_st.k = 16.2; % [W/mK]
const_st.cp = 500; % [J/kgK]
const_st.rho = 8000; % [kg/m^3]
const_st.alpha = const_st.k/(const_st.rho*const_st.cp); % [m^2/s]
const_st.L = x8+0.0508; % [m]

const_br26.H_exp = 114; % [C/m]
const_br26.T0 = 12; % [C]
const_br26.k = 115; % [W/mK]
const_br26.cp = 380; % [J/kgK]
const_br26.rho = 8500; % [kg/m^3]
const_br26.alpha = const_br26.k/(const_br26.rho*const_br26.cp); % [m^2/s]
const_br26.L = x8+0.0508; % [m]
 
const_br29.H_exp = 139.8; % [C/m]
const_br29.T0 = 12; % [C]
const_br29.k = 115; % [W/mK]
const_br29.cp = 380; % [J/kgK]
const_br29.rho = 8500; % [kg/m^3]
const_br29.alpha = const_br29.k/(const_br29.rho*const_br29.cp); % [m^2/s]
const_br29.L = x8+0.0508; % [m]

const_al26.H_exp = 54.3; % [C/m]
const_al26.T0 = 12; % [C]
const_al26.k = 130; % [W/mK]
const_al26.cp = 960; % [J/kgK]
const_al26.rho = 2810; % [kg/m^3]
const_al26.alpha = const_al26.k/(const_al26.rho*const_al26.cp); % [m^2/s]
const_al26.L = x8+0.0508; % [m]

const_al28.H_exp = 66.8; % [C/m]
const_al28.T0 = 12; % [C]
const_al28.k = 130; % [W/mK]
const_al28.cp = 960; % [J/kgK]
const_al28.rho = 2810; % [kg/m^3]
const_al28.alpha = const_al28.k/(const_al28.rho*const_al28.cp); % [m^2/s]
const_al28.L = x8+0.0508; % [m]

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


for j = 1:length(alpha_st)
    const_st.alpha = alpha_st(j);
    for i = 1:length(t_st)
    u_st(i,:) = heatdistribution(t_st(i),x,const_st);
    end
    error_st(j) = mean(rmse(u_st, st_notime));
end
[~,idx] = min(error_st);
const_st.alpha = alpha_st(idx);

for j = 1:length(alpha_al26)
    const_al26.alpha = alpha_al26(j);
    for i = 1:length(t_al26)
    u_al26(i,:) = heatdistribution(t_al26(i),x,const_al26);
    end
    error_al26(j) = mean(rmse(u_al26, al26_notime));
end
[~,idx] = min(error_al26);
const_al26.alpha = alpha_al26(idx);
    
for j = 1:length(alpha_al28)
    const_al28.alpha = alpha_al28(j);
    for i = 1:length(t_al28)
    u_al28(i,:) = heatdistribution(t_al28(i),x,const_al28);
    end
    error_al28(j) = mean(rmse(u_al28, al28_notime));
end
[~,idx] = min(error_al28);
const_al28.alpha = alpha_al28(idx);

for j = 1:length(alpha_br26)
    const_br26.alpha = alpha_br26(j);
    for i = 1:length(t_br26)
    u_br26(i,:) = heatdistribution(t_br26(i),x,const_br26);
    end
    error_br26(j) = mean(rmse(u_br26, br26_notime));
end
[~,idx] = min(error_br26);
const_br26.alpha = alpha_br26(idx);

for j = 1:length(alpha_br29)
    const_br29.alpha = alpha_br29(j);
    for i = 1:length(t_br29)
    u_br29(i,:) = heatdistribution(t_br29(i),x,const_br29);
    end
    error_br29(j) = mean(rmse(u_br29, br29_notime));
end
[~,idx] = min(error_br29);
const_br29.alpha = alpha_br29(idx);

figure()
hold on
grid minor
plot(alpha_st,error_st,'k')
scatter(const_st.alpha,min(error_st),'filled','r')
xlabel('\alpha [m^2/s]')
ylabel("RMSE")
title('RMSE vs \alpha: Steel - 21V,192mA')

figure()
hold on
grid minor
plot(alpha_al26,error_al26,'k')
scatter(const_al26.alpha,min(error_al26),'filled','r')
xlabel('\alpha [m^2/s]')
ylabel("RMSE")
title('RMSE vs \alpha: Aluminum - 26V,250mA')

figure()
hold on
grid minor
plot(alpha_al28,error_al28,'k')
scatter(const_al28.alpha,min(error_al28),'filled','r')
xlabel('\alpha [m^2/s]')
ylabel("RMSE")
title('RMSE vs \alpha: Aluminum - 28V,269mA')

figure()
hold on
grid minor
plot(alpha_br26,error_br26,'k')
scatter(const_br26.alpha,min(error_br26),'filled','r')
xlabel('\alpha [m^2/s]')
ylabel("RMSE")
title('RMSE vs \alpha: Brass - 26V,245mA')

figure()
hold on
grid minor
plot(alpha_br29,error_br29,'k')
scatter(const_br29.alpha,min(error_br29),'filled','r')
xlabel('\alpha [m^2/s]')
ylabel("RMSE")
title('RMSE vs \alpha: Brass - 29V,273mA')

for i = 1:length(t_st)
    u_st(i,:) = heatdistribution(t_st(i),x,const_st);
end
for i = 1:length(t_al26)
    u_al26(i,:) = heatdistribution(t_al26(i),x,const_al26);
end
for i = 1:length(t_al28)
    u_al28(i,:) = heatdistribution(t_al28(i),x,const_al28);
end
for i = 1:length(t_br26)
    u_br26(i,:) = heatdistribution(t_br26(i),x,const_br26);
end
for i = 1:length(t_br29)
    u_br29(i,:) = heatdistribution(t_br29(i),x,const_br29);
end

figure()
hold on
for i = 1:8
plot(al26_t,al26_c{i},'k')
plot(t_al26,u_al26(:,i),'r')
end
grid minor
xlabel("Time (s)")
ylabel('Temperature (^\circC)')
legend('Experimental','Model','','','','','','','','','','','','','','','Location','northwest')
title("Model III: Aluminum - 26V,250mA")

figure()
hold on
for i = 1:8
plot(al28_t,al28_c{i},'k')
plot(t_al28,u_al28(:,i),'r')
end
grid minor
xlabel("Time (s)")
ylabel('Temperature (^\circC)')
legend('Experimental','Model','','','','','','','','','','','','','','','Location','northwest')
title("Model III: Aluminum - 28V,269mA")

figure()
hold on
for i = 1:8
plot(br26_t,br26_c{i},'k')
plot(t_br26,u_br26(:,i),'r')
end
grid minor
xlabel("Time (s)")
ylabel('Temperature (^\circC)')
legend('Experimental','Model','','','','','','','','','','','','','','','Location','northwest')
title("Model III: Brass - 26V,245mA")

figure()
hold on
for i = 1:8
plot(br29_t,br29_c{i},'k')
plot(t_br29,u_br29(:,i),'r')
end
grid minor
xlabel("Time (s)")
ylabel('Temperature (^\circC)')
legend('Experimental','Model','','','','','','','','','','','','','','','Location','northwest')
title("Model III: Brass - 29V,273mA")

figure()
hold on
for i = 1:8
plot(st_t,st_c{i},'k')
plot(t_st,u_st(:,i),'r')
end
grid minor
xlabel("Time (s)")
ylabel('Temperature (^\circC)')
legend('Experimental','Model','','','','','','','','','','','','','','','Location','northwest')
title("Model III: Steel - 21V,192mA")


function u = heatdistribution(t,x,const)
    sum = [0,0,0,0,0,0,0,0];
    for n = 1:10
        b_n = (8*const.H_exp*const.L*(-1)^n)/((2*n-1)^2*pi^2);
        lambda_n = pi*(2*n-1)/(2*const.L);
        sum = sum + b_n.*sin(lambda_n.*x).*exp(-1*lambda_n^2*const.alpha*t);
    end
    u = const.T0+const.H_exp.*x + sum;
end
