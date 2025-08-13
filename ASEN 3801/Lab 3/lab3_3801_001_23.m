close all; clear; clc;

% Storing Data
data{1,1} = readmatrix('2024_10_04_001_RWHEEL_1');
data{2,1} = readmatrix('2024_10_04_001_RWHEEL_i01_f01_t60s');
data{3,1} = readmatrix('2024_10_04_001_RWHEEL_I01'); %f=1
data{4,1} = readmatrix('2024_10_04_001_RWHEEL_I1_f01_t60s');
data{5,1} = readmatrix('2024_10_04_001_RWHEEL_I1_f08_60s');

%Preallocation
time = cell(5,1);
controlGyro = cell(5,1);
bus = cell(5,1);
coefficients = cell(5,1);
xFit = cell(5,1);
yFit = cell(5,1);
b = cell(5,1);
calibrated = cell(5,1);
CalibratedError = cell(5,1);
rawError = cell(5,1);
posTrue = cell(5,1);
gyroPos = cell(5,1);
posError = cell(5,1);


% Title string vector for the later plots
titles = ["Test Example","Current = 0.1 A, Frequency = 0.1 Hz","Current = 0.1 A, Frequency = 1 Hz","Current = 1 A, Frequency = 0.1 Hz","Current = 1 A, Frequency = 0.8 Hz"];
for i = 4

    % Seperate data into columns, do some conversions
    time{i} = data{i}(2 :end,1);
    t = ones(length(time{i}),1);
    time{i} = time{i}(:,1)-time{i}(1,1).*t;
    controlGyro{i} = data{i}(2: end,2); %rad/s
    bus{i} = data{i}(2: end,3).*pi/30; %rad/s
    
    %% a) Time history of measurements

% Subplot Graphs
%     figure()
%     subplot(2,1,1)
%     plot(time{i},bus{i},'Linewidth',1)
%     title('Bus angular velocity (Encoder Measurement)')
%     grid minor;
%     xlabel('Time (s)');
%     ylabel('Angular Velocity (rad/s)')
% 
%     subplot(2,1,2)
%     plot(time{i},controlGyro{i},'Linewidth',1)
%     title('Control Gyro Angular Velocity');
%     grid minor;
%     xlabel('Time (s)');
%     ylabel('Angular Velocity (rad/s)')
%     sgtitle(sprintf(titles(i)));

% Plot with the graphs overlaid
    figure()
    hold on;
    plot(time{i},bus{i},'linewidth',1);
    plot(time{i},-controlGyro{i},'LineWidth',1)
    grid on;
    xlabel("Time (s)");
    ylabel("Component Angular Velocity (rad/s)")
    legend("Bus","Reaction Wheel");
    title(sprintf('Comparison of Reaction Wheel and Bus Angular Velocity for %s',titles(i)))

    %% b) Time history of angular rate measurement error
% Error from encoder
    %     figure();
%     subplot(2,1,1);
%     hold on;
%     scatter(bus{i},controlGyro{i},'Marker','.');
%     grid minor;
%     title('Encoder Measurement vs Gyro measurement, Uncalibrated');
%     xlabel("Encoder Rate (rad/s)");
%     ylabel("Gyro Rate Measurement (rad/s)");
    
    % Calculating the line of best fit of the correlation
    coefficients{i} = polyfit(bus{i}, controlGyro{i}, 1);
    % Create a new x axis with exactly 1000 points 
    xFit{i} = linspace(min(bus{i}), max(bus{i}), 1000);
    % Get the estimated yFit value for each of those 1000 new x locations.
    yFit{i} = polyval(coefficients{i} , xFit{i});

    % Correlation best fit plot
%     plot(xFit{i},yFit{i},'r--','Linewidth',1.5);
%     yline(coefficients{i}(2),'linewidth',1)
%     xlabel('Encoder Measurement (rad/s)');
%     ylabel('Calibrated Gyro Measurement (Rad/s)');
%     legend("Data",'Correlation','Bias b')

    % Calibrating data using coeffients from best fit
    b{i} = coefficients{i}(2)*ones(length(time{i}),1);
    calibrated{i} = 1/coefficients{i}(1)*(controlGyro{i}-b{i});

    % Correlation of calibrated data
%     subplot(2,1,2);
%     scatter(bus{i},calibrated{i},'Marker','.');
%     grid minor;
%     xlabel('Encoder Rate Measurement (rad/s)')
%     ylabel('Calibrated Gyro Rate (rad/s)')
%     title('Calibrated Data')
% 
%     sgtitle(sprintf('Encoder Gyro Measurement Correlation for %s',titles(i)))


%     figure()
%     hold on;
%     plot(time{i},bus{i},'linewidth',1);
%     plot(time{i},calibrated{i},'LineWidth',1)
%     grid on;
%     xlabel("Time (s)");
%     ylabel("Component Angular Velocity (rad/s)")
%     legend("Bus","Reaction Wheel");
%     title(sprintf('Comparison of Calibrated Reaction Wheel and Bus Angular Velocity for %s',titles(i)))

    %% c) Time history of angular position error

    % Calculating the Error in Measurement
    CalibratedError{i} = abs(bus{i})-abs(calibrated{i});
    rawError{i} = abs(controlGyro{i}) - abs(bus{i});

    meanErrorRaw(i) = mean(rawError{i});
    meanErrorCalibrated(i) = mean(CalibratedError{i});

    figure();
    hold on
    subplot(2,1,1);
    scatter(time{i}, rawError{i},'marker','.')
    title('Measurement Error vs. Time')
    yline(meanErrorRaw(i),'Linewidth',1)
    xlabel('Time (s)');
    ylabel('Error (rad/s)');
    grid minor;
    legend('','Mean Error','location','southwest')
    hold off;

    subplot(2,1,2);
    hold on;
    scatter(time{i}, CalibratedError{i},'marker','.')
    title('Calibrated Data Error vs. Time')
    yline(meanErrorCalibrated(i),'Linewidth',1)
    xlabel('Time (s)');
    ylabel('Error (rad/s)');
    grid minor;
    legend('','Mean Error','location','southwest')
    sgtitle(sprintf('Time History of Angular Rate Error (measured vs. calibrated) for %s',titles(i)))
    hold off;

    %% d) Time history of angular position error
    deltaT = time{i}(2)-time{i}(1);
    intBus = deltaT*bus{i};
    intGyro = deltaT*controlGyro{i};

    for j = 1:length(time{i})
        posTrue{i}(j,1) = sum(intBus([1,j]));
        gyroPos{i}(j,1) = sum(intGyro([1,j]));

    end

    figure()
    subplot(3,1,1);
    plot(time{i},posTrue{i},'Linewidth',1);
    grid minor;
    xlabel('Time (s)');
    ylabel('Angular Position (rad)');
    title('Angular Position from Encoder');

    subplot(3,1,2);
    plot(time{i},gyroPos{i},'Linewidth',1);
    grid minor;
    xlabel('Time (s)');
    ylabel('Angular Position (rad)');
    title('Angular Position from Gyro');
    
    sgtitle(sprintf('Time History of Angular Position Error %s',titles(i)));

    posError{i} = gyroPos{i} + posTrue{i};
    subplot(3,1,3);
    plot(time{i},posError{i},'Linewidth',1)
    grid minor;
    xlabel('Time (s)');
    ylabel ('Angular Error (rad)')
    title('Position Error')

    figure()
    subplot(2,1,1);
    plot(time{i},posTrue{i},'Linewidth',1);
    grid minor;
    xlabel('Time (s)');
    ylabel('Angular Position (rad)');
    title('Angular Position from Encoder');

    subplot(2,1,2);
    plot(time{i},gyroPos{i},'Linewidth',1);
    grid minor;
    xlabel('Time (s)');
    ylabel('Angular Position (rad)');
    title('Angular Position from Gyro');
    sgtitle()
    sgtitle(sprintf('Time History of Angular Position %s',titles(i)));


    figure()
    plot(time{i},posError{i},'Linewidth',1)
    grid minor;
    xlabel('Time (s)');
    ylabel ('Angular Error (rad)')
    title('Position Error')

    %% e) Angular Position Error vs encoder rate

    figure()
    scatter(bus{i},posError{i},'marker','.');
    grid minor;
    title(sprintf('Angular Position Error vs Encoder Rate for %s',titles(i)));
    ylabel('Angular Position Error (rad/s)')
    xlabel('Encoder Rate Measurement (rad/s)');


end

b_mean = b{1,1}(1,1)+b{2,1}(1,1)+b{3,1}(1,1)+b{4,1}(1,1)+b{5,1}(1,1);

bias = zeros(5,1);
sensitivity = zeros(5,1);
for i = 1:5

    bias(i) = coefficients{i}(2);
    sensitivity(i) = coefficients{i}(1);

end

[devBias,mBias] = std(bias)
[devSens,mSens] = std(sensitivity)



