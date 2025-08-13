clc;close all;clear;

 

data{1,1} = load("2024_10_04_001_RWHEEL_T5t3");
data{2,1} = load("2024_10_04_001_RWHEEL_T5t5");
data{3,1} = load("2024_10_04_001_RWHEEL_T7t3");
data{4,1} = load("2024_10_04_001_RWHEEL_T7t5");
data{5,1} = load("2024_10_04_001_RWHEEL_T7t7");

 

tau = 33.5/1000;     %From data sheet, Nm/A, using motor #251601
rpm = pi/30;    %Conversion for RPM to rad/s
T = [5,5,7,7,7] / 1000;    %Torque values from file names in Nm
t = [3,5,3,5,7];    %Time values from files
Aero_torque = 10^-4; %Torque for wheel to resist
wheel_max = rpm*4000;   %max spin rate in rad/s
tstart = 110;   %This is the approximate index for where the wheel begins accelerating
r_wheel = 42.8 / 1000; %Wheel radius in meters from data sheet

 

%Preallocating cells for for loop speed

time = cell(5,1);
torque = cell(5,1);
current = cell(5,1);
omega = cell(5,1);
coeffs = cell(5,1);
yfit = cell(5,1);
MOI = cell(5,1);

 

for i=1:5

    current{i,1} = data{i}(:,4);    %Current used for torque
    torque{i,1} = tau*current{i,1}; %Converting current to torque [Nm]
    omega{i,1} = rpm*data{i}(:,3);  %Angular velocity, converted to rad/s

 
    time{i,1} = data{i}(:,1)/1000;  %Time for each file in seconds
    [M(i),I(i)] = max(omega{i});
    time{i} = time{i}(tstart:I(i));
    omega{i} = omega{i}(tstart:I(i));
    torque{i} = torque{i}(tstart:I(i));

 

    %Line of best fit
    coeffs{i} = polyfit(time{i},omega{i},1);
    slope(i) = coeffs{i}(1);   %Average acceleration, rate of change of best fit
    slope(i) = slope(i)*r_wheel;    %Converting angular acceleration to linear acceleration
    yfit{i} = polyval(coeffs{i},time{i});

    %Plotting for each file
    figure(i)
    hold on;
    grid on;
    plot(time{i},omega{i})
    plot(time{i},yfit{i})
    xlabel('Time [s]')
    ylabel('Angular Velocity [rad/s]')
    legend('Angular Velocity','Line of Best Fit','location','best')
    title(sprintf('Angular Velocity vs Time for Torque of %g mNm and Time of %g Seconds',T(i),t(i)))

 
    MOI{i} = mean(torque{i}) / slope(i);  %Calculating moment of inertia

end

 

%Calculating average moment of inertia

sum = 0;

for i = 1:5

    sum = sum + MOI{i};
    numb(i) = MOI{i};

end

 

Avg_Inertia = sum/5;
Standard_Dev_Inertia = std(numb);


fprintf("The mean moment of inertia is %g kgm^2 with a standard deviation of %g. \n",Avg_Inertia,Standard_Dev_Inertia);

 

%% Finding the time before wheel rpm exceeds 4000 with T=10^(-4) Nm

 

%Using third case as basis

 

accel = slope(3);

omega = 4000*rpm;