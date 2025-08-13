clear;
clc;
close all; 
data=readmatrix('Static Test Stand Calibration Case 2.xlsx'); %read the data from the file called Static Test Stand Calibration Case 2.xlsx
weight=data(:,1);% the first column of the matrix data
f0_offset=data(:,2);% the second column of the matrix data 
f1_offset=data(:,3); % the third column of the matrix data 
f0_mv=data(:,4); %the fourth column of the matrix data 
f1_mv=data(:,5);% the fifth column of the matrix data 
range=(1:1:length(data));
%calculate forces load weight for f1 and f0
f0_ef=f0_mv-f0_offset;
f1_ef=f1_mv-f1_offset;
Sort1=sort(f0_ef);
Sort2=sort(f1_ef);

[P,s]= polyfit(range,Sort1,1); 
%fitting parameters 
a0=P(1); % slope
a1=P(2); %intercept
slope_fit=a1+a0.*range; % slope of best fit line 

%plot with the original linear scale
figure();
hold on;
scatter(range,Sort1);
plot(range,slope_fit,'b-')
%labels 
xlabel('Range');
ylabel('Force of F0');
title('Best Fit Line of F0');
text(23,5,string);
sprintf('a0= %2.4f, a1= %2.4f',a0,a1);%Add text to the plot
grid on;%diplay grid lines
[R,T]= polyfit(range,Sort2,1); 
%fitting parameters 
a2=R(1); % slope
a3=R(2); %intercept
slope_fit2=a3+a2.*range;

%plot with the original linear scale
figure();
hold on;
scatter(range,Sort2);
plot(range,slope_fit2,'b-');
%labels 
xlabel('Range');
ylabel('Force of F1');
title('Best Fit Line of F1');
text(23,5,string);
sprintf('a0= %2.4f, a1= %2.4f',a2,a3);%Add text to the plot
grid on;%diplay grid lines






