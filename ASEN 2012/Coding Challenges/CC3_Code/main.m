% Author(s): Andrew Patella
% Assignment title: Coding Challenge 3
% Purpose:
% Creation date: 10/16/2023
% Revisions: N/A

clear; clc; close all;

%Read in data
data = readmatrix("supercomputer-power-flops.csv");

%Storing data in matrices
time = data(:,1);
time = time-1990; %years
flops = data(:,2); %GFLOPS

%Creating an arbitrary time array that extends outside the time data bounds
time_long = -5:1:40;

%P value for exponential regression
p= 2;

%Getting polynomial coefficients from polyfit, first degree.
[line,S] = polyfit(time,log(flops),1);

%Changing linear polynomials to exponential polynomials
a_1 = exp(line(2));
a_0= line(1)/log(p);

%Function that produces line of best fit, working on whatever time interval
fit =@(t) a_1*p.^(a_0*t);

%Creating a line of best fit for data time and extended time
fit_1 = fit(time);
fit_long = fit(time_long);

%string storing a_0 and a_1 values
string = sprintf('a_0 = %2.4f, a_1 = %2.4f',line(1),line(2));

%Using polyval to find delta and extrapolated fit data
[extrap_fit,delta] = polyval(line,time_long,S);

%Creating arrays of error bars for calculation
top_error=exp(extrap_fit+2*delta);
bottom_error = exp(extrap_fit-2*delta);

%plotting the results
figure(1);
hold on;
scatter(time,flops);
set(gca,'YScale','log');
plot(time_long,fit_long,'linewidth',1,'color','blue');
text(21,fit(31.5),string,'FontAngle','italic');
xlabel('Years since 1990');
ylabel('Number of Flops');
grid on;
plot(time_long,exp(extrap_fit+2*delta),'linewidth',.75,'color','red','LineStyle','--');
plot(time_long,exp(extrap_fit-2*delta),'linewidth',.75,'color','red','LineStyle','--');
title('Supercomputer GFLOPs by year');
legend('Data','line of best fit','\pm 2\sigma_q','Location','southeast');

hold off;

%Calculating and printing differences in error bars
gflops_1990 = fit_long(6);
error_beginning = (top_error(6))-(bottom_error(6));
error_2025 = (top_error(41))-(bottom_error(41));
fprintf('The estimated number GFLOPS in 1990 is %2.4f\n',gflops_1990);
fprintf('The difference in the error bars (95%% confidence) is %2.1f GFLOPS in 1990 \n',error_beginning);
fprintf('The difference in the error bars (95%% confidence) is %2.2g GFLOPS in 2025 \n',error_2025);


%% Comments 

% The prediction I would trust the most is the one for 1985. This one is
% only 5 years away from the initial known value. Therefore there is less error. This can also be seen in how the the 95% error interval is smaller at this time.  
% In 2025, the years are further away so the error bounds are growing away
% from the best fit line. There is less certainty in extrapolation the
% further away from the data. 

