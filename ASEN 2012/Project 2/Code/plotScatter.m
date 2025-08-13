function [trajectory_max,ideal] = plotScatter(test_parameter,state,parameter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function is designed to create scatterplots to demonstrate trends
% when initial conditions of the rocket are changed. This allows for
% optimization even if a direct function cannot be calculated.

% This function also creates a line of best fit for the data. This allows
% it to return a more reasonable max value for the optimal value without
% being effected by any outliers in the data. It creates a 6th degree
% polynomial for best fit. 

% INPUTS:

% test_parameter is the array of different values that are being checked

% state is the cell with corresponding state of the rocket from ode45 for
% each of the values of the changed parameter

% parameter is a string for the labels and title of the scatter plot

% OUTPUTS: 

% trajectory_max is a value how far the rocket could go with the ideal
% condition

% ideal is the ideal paramater value for maximum distance.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %This is an array of the max distance value for each row of the state
    %cell
    %Preallocation
    maxD = zeros(length(test_parameter),1);
    maxi = zeros(length(test_parameter),1);
    maxH = zeros(length(test_parameter),1);
    maxIh = zeros(length(test_parameter),1);

    %Calculation
    for i = 1: length(test_parameter)
        %The max and when of each matrix in the cell
        [maxD(i),maxi(i)] = max(state{i}(:,1));
        [maxH(i),maxIh(i)] = max(state{i}(:,3));
    end

    %Calculating the 6th degree polynomial of best fit in order to
    %reasonabally ignore any outliers
    [P,S,Mu] = polyfit(test_parameter,maxD,6);
    [Ph,Sh,Muh] = polyfit(test_parameter,maxH,6);
    [Y,delta] = polyval(P,test_parameter,S,Mu);
    [Yh,deltah] = polyval(Ph,test_parameter,Sh,Muh);

    %Assuming the best fit is a reasonable fit, the max of this line is
    %approximately the max of distance that can be traveled with the
    %changed parameter
    [trajectory_max(1),index] = max(Y);
    [trajectory_max(2),indexh] = max(Yh);

    %Calculating the ideal value from the index in the cell that
    %corresponds to maximum distance
    ideal(1) = test_parameter(index);
    ideal(2) = test_parameter(indexh);


    %Creating a scatter plot
    figure()
    hold on;
    scatter(test_parameter,maxD);
    grid on;
    xlabel(sprintf(parameter));
    ylabel('Corresponding max distance');
    title(sprintf('Max distance due to changing %s',parameter));
    plot(test_parameter,Y);
    plot(test_parameter,Y+2*delta);
    plot(test_parameter,Y-2*delta);

    figure()
    hold on;
    scatter(test_parameter,maxH);
    grid on;
    xlabel(sprintf(parameter));
    ylabel('Corresponding max Height');
    title(sprintf('Max Height due to changing %s',parameter));
    plot(test_parameter,Yh);
    plot(test_parameter,Yh+2*deltah);
    plot(test_parameter,Yh-2*deltah);
    



end