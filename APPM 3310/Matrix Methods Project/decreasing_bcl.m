%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script focuses on recreating multiple pulses for different
% heartbeats and analyzing the eigenvalues of that system. The heartbeats
% are synthesized by random numbers and are inherently biased but should be
% releatively useful.

% THIS ONE IS DIFFERENT BECAUSE IT VARIES CYCLE LENGTH (BCL)

% LIST OF ASSUMPTIONS THAT MAY BE (PROBABLY ARE) FALSE

% - Steady state value can be approximated by the mean
% - Every beat is the same length
% - The data we have is reasonably accurate

% Change first and last to make sure that the APD values are matching- also
% divide by ten. Change the number of heartbeats to make sure you are
% matching their assumptions. Currently at n=7 beats (>6) and 6 APD values
% from 40 to 90 (first = 4, last = 9). 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clc; clear;
n = 4; % Number of heartbeats (this will graphically add one more no se)

% APD^[] values that you wish to measure (x/10)
first = 4;
last = 9;

% Number of Dominant Compnents
m = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Read in data from last script which conditioned the data to be within
%reasonable values. From "matrix_project.m"
tmvm = readmatrix("tmv.csv");
time = readmatrix("time.csv");
APD_ss = readmatrix("APD_ss");

time = time/1;


values = first:1:last;
tmv_tot = tmvm;
t_tot = time;
duration_beat = time(end);
endIndex = 1;
magnitude = abs(max(tmvm)-min(tmvm));
num_vals = last-first;

%Prealocating
Z = zeros(n,num_vals);
peakV = zeros(n,1);
peaki = zeros(n,1);
riseV = cell(n,1);%
riseT = cell(n,1);
fallT = cell(n,1);
fallV = cell(n,1);
t2 = cell(last,1);
tmv = cell(1,n);
v_apd = cell(1,last);
v_diff = zeros(n,1);
APDsc = cell(1,last);
APDs = zeros(num_vals,1);
t_stim = zeros(1,n);
APDv = cell(1,last);
Y1 = zeros(n,m);
Y2 = zeros(n,m);
BCL = zeros(n-1,1);




%% Creating the random heartbeats

for i = 1:n
    for j = 1:size(tmvm)
        tmv{i}(j,1) = tmvm(j)+0.05*abs(magnitude+tmvm(j))*randn(1,1);
        t{i}(j,1) = time(j)+ 10*randn(1,1)^2;
        %t{i}(j,1) = time(j) - time(j)/1.5;
    end

    time = time + t{i}(end,1) - t{i}(1,1);
    tmv_tot = cat(1,tmv_tot,tmv{i});
    t_tot = cat(1,t_tot,time);
end

figure()
hold on;
plot(t_tot,tmv_tot,"Linewidth",1);
xlabel("Time [ms]");
ylabel("Transmembrane Voltage [mV]");
title(sprintf("Transmembrane Voltage for %g beats",n+1));
grid minor;



%% Analyzing TMV data for APD values


% Outer loop iterates through each heartbeat
for j = 1: n

    %Splitting the individual beats into rise and fall
    [peakV(j),peaki] = max(tmv{j});

    riseV{j} = tmv{j}(endIndex:peaki);
    riseT = time(endIndex:peaki);

    fallV{j} = tmv{j}(peaki+9:end);
    fallT = time(peaki+9:end);
 
    endIndex = endIndex+ size(tmvm,1);

    v_diff(j) = abs(max(tmv{j})-min(tmv{j}));


    % Iterating through the different values of APD desired
    for i = 1:length(values)
        v_apd{values(i)} = peakV(j)-values(i)/10*v_diff(j);
    
        t2{values(i)} = interp1(fallV{j},fallT,v_apd{values(i)});

        t_stim(j) = t_tot(endIndex+1);
        
        APDv{values(i)} = t2{values(i)}-t_stim;
    end
end


%% Forming the matrix, SVD

y_star = APDs;


for i = 1:n
    for j = first:last
        Z(i,j-(first-1)) = APDv{j}(i)-APD_ss(j-(first-1)); % i is beat number, j is APD value
    end
end

Z = flip(Z,2);


% Singular Value Decomposition of Z
[U,S,V] = svd(Z,0);  % MAYBE ACTUALLY DO THIS? Z = USV'


% zb should be mxN, N=n=#beats
Zb = Z(:,1:m);

%Calculating Y1 and Y2
for i = 1:n-1
    Y1(i,:) = Zb(i,:)';
end

for i = 2:n
    Y2(i,:) = Zb(i,:)';
end

Y1 = Y1';
Y2 = Y2';

%% Calculating D
D = (Y2*Y2')*inv(Y1*Y1');

% Calculating the eigenvalues of the system
eigenvalue = eig(D)

for i = 1:size(t_stim,2)-1

    BCL(i) = t_stim(i+1)-t_stim(i);
end

BCL = mean(BCL)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



