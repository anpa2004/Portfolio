% Author(s): Andrew Patella
% Assignment title: Coding Challenge 2
% Purpose: Compare Monte Carlo and General Method models for error propagation
% Creation date: 10/2/2023
% Revisions: N/A

clear;
clc;
close all;

%% Read in and Clean the Data Set

%Storing data from file in array
table_data = readmatrix('Static_Thrust_Data.xlsx');

%Cleaning out all thrust data where thrust was less than 100 N or nan
Thrust1 = table_data(:,2);
Thrust1(Thrust1<1000) = [];
Thrust1(isnan(Thrust1)) = [];

Thrust2 = table_data(:,3);
Thrust2(Thrust2<1000) = [];
Thrust2(isnan(Thrust2)) = [];

%cleaning up the time column to be the length of the thrust columns with
%the same time incriment as the thrust columns
time = 0:.05:(length(Thrust1)*.05)-.05;
time = transpose(time);

%Combining the thrust columns into one array. 
thrust(:,1) = Thrust1(:,1);
thrust(:,2) = Thrust2(:,1);

thrustAvg = sum(thrust(:))/length(thrust(:));
thrustStd = std(thrust(:));

%% Calculated Weighted average and Uncertainty
%Inputting the thrust matrix to the matrix_std function which outputs
%average and the uncertainty
[avgThrust,uncertainty] = matrix_std(thrust);

%Calculating the mean and standard deviation for the thrust
mThrust = mean(thrust(:));
sdThrust = std(thrust(:));

%Plotting the thrust data and mean and standard deviation
figure(1);
hold on;
plot(time,thrust(:,1));
plot(time,thrust(:,2));
yline(mThrust,'linewidth',1);
yline(mThrust+sdThrust,'--','linewidth',1);
yline(mThrust-sdThrust,'--','linewidth',1);
grid on;
legend('Sensor 1 Thrust data','Sensor 2 Thrust data','Mean Thrust','Mean thrust \pm \sigma_x');
xlabel('Time (s)');
ylabel('Thrust (N)');
title('Thrust Profile vs. Time');
hold off;

%% The General Method

% Declaring constants
g_0 = 9.81; %m/s^2

%I_sp and its uncertainty
I_sp = 459;  %s
dI_sp = 11;
%M_s and its uncertainty
m_s = 13050; %kg
dm_s = 60; %kg
%M_p and its uncertainty
m_p = 71800; %kg
dm_p = 300; %kg
%M_0 and its uncertainty
m_0 = m_s + m_p;
dm_0 = dm_s+dm_p;
%m_f and its uncertainty
m_f = m_s;
dm_f=dm_s;

%Partial derivatives of the change in velocity formula
dvdI = g_0*log(m_0/m_f)*dI_sp;
dvdm_0 = I_sp*g_0*dm_0/m_0;
dvdm_f =-I_sp*g_0*dm_f/m_f;

%Delta v formula
dv = sqrt(dvdI^2+dvdm_f^2+dvdm_0^2);
%V formula
v = I_sp*g_0*log(m_0/m_f);

%Outputting the change in velocity and its uncertainty
fprintf('Change in velocity is %4.0f \x00B1%4.0f m/s',v,dv);

%% Monte Carlo Simulation

%Declaring number of random values for monte carlo
N = 20000;

% Declaring variables, creating a matrix of N possible values
I = 459 + 11.*randn(1,N);
MCm_p = 71800+300.*randn(1,N);
MCm_f = 13050+60.*randn(1,N);
MCm_0 = MCm_p+MCm_f;

%calculating the change in velocity 
vMC = g_0*I*log(MCm_0/MCm_f);


%Calculating mean and std of the velocity matrix
MCmean =mean(vMC);
MCstdev = std(vMC);


%% Compare Monte Carlo and General Method
[PDF, x] = create_normal_distribution(v,dv,min(vMC),max(vMC),N);

% Plotting the results
figure(2);
hold on;
histogram(vMC,500);
%Mean from monte carlo, with +-1 standard deviation
xline(MCmean,'r--','linewidth',1);
xline(MCmean-MCstdev,'--');
xline(MCmean+MCstdev,'--');

%Labels and title
xlabel('Change in Velocity (m/s)');
ylabel('Occurances');
title('Change in velocity for Monte Carlo and General Method error propagation methods');

%Right y axis so the pdf is visible
yyaxis right;
plot(x,PDF,'linewidth',1,'color','r');
grid on;


legend('Data from Monte Carlo','Mean \Delta v (Monte Carlo)','\Delta v \pm \sigma_v (Monte Carlo)','','PDF (General Method)');




function [pdf,x] = create_normal_distribution(xMean,xStd,xMin,xMax,N)
%create_normal_distribution Create normal distribution vector 
%   calculates a normal distribution from the minimum to maximum with N
%   components 
x=linspace(xMin,xMax,N);
pdf = 1/(sqrt(2*pi*xStd.^2)) .* exp (- (x-xMean).^2 ./ (2*xStd.^2));
end

%% Questions

% The Monte Carlo simulation and general methods are very similar with their 
% estimates. They both provide slightly different valus but when plotted 
% with the histogram for the Monte Carlo, the solutions are very similar. If 
% considering the uncertainty with 2*sigma, the uncdrtainty will increase
% and there will be more difference between the models. The accepted values
% for delta_v increases so the Monte Carlo would have a greater range than
% the general. Calculating the average and standard deviations for the
% thrust profiles and using weighted averages yields different values. They
% are similar but not exact. The weighted average takes into account the
% uncertainty of specific measurements. The average is just assuming the
% same certainty for all measurements. 
