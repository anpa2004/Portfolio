% Author(s): Andrew Patella, Austin Vera
% Assignment title: Project 2
% Purpose: Model how changing initial parameters changes the trajectory of
% a bottle rocket
% Creation date: 11/06/2023
% Revisions: N/A

close all; clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script calculates the trajectories for different initial
% parameters.

% Initially, it uses 8 points for plotting trajectories. This is because
% too many becomes cluttered on the plots. 

% The second part calculates trajectories for n parameter values between
% initial and final values for a scatter plot. More points means more
% certainty in the prediction of optimal parameter value because it
% has more resolution so the approsimation is more accurate. 

% The third part of this function is a process of refining initial
% conditions to get the trajectory as close as possible to 80 m. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Constants and initial conditions

% Storing constants in a struct using the getConst() function
const = getConst();

%Integration time vector
int_time = [0 5];

%% Change in performance due to change in parameters

% Arrays of 8 values of each parameter to check
% Only 8 were chosen initially just to provide helpful plots that aren't
% too crowded. 
theta_test = [25,30,35,40,45,50,55,60]; %Degrees
theta_test = theta_test*pi/180; %Radians

pressure_test = [45,50,55,60,65,70,75,80]; %psi
pressure_test = pressure_test*6894.76; %Pa

V_test = [0.001,0.1,0.2,0.3,.4,.5,.6,.65]; %Percent
V_test = const.V_b * V_test; %Volume

cd_test = [0.3,0.35,0.38,0.4,0.41,0.45,0.48,0.5];

%Using local varying_parameters function to create a cell of the states
%from state matrix func for each to see differences due to a new parameter
[~,state_theta] = varying_parameters("theta_i",theta_test,int_time);
[~,state_P] = varying_parameters("p_0",pressure_test,int_time);
[~,state_c] = varying_parameters("c_D",cd_test,int_time);
[~,stateV] = varying_parameters("V_0w",V_test,int_time);

%Plotting the trajectories with corresponding changed values
plotFun(state_theta,'theta_i');
plotFun(state_P,'p_0');
plotFun(state_c,'c_D');
plotFun(stateV,'V_0w');


%% Scatter Plots for optimization

%Max and min initial pressures (MUST BE IN PSI)
pmax = 80; 
pmin = 45;

%Max and min coefficients of drag (unitless)
cmax = 0.3;
cmin = 0.5;

%Max and min launch angle (MUST BE IN DEGREES)
thetamin = 0;
thetamax = 90;

%Max and min percents of water
Vmin = 0.001;
Vmax = 0.6;

% Number of values, times ode45 will run with different parameters
npts = 500;

% Long matrices between max and min values for more scatter plot (Accuracy)
p_long = linspace(pmin,pmax,npts);
p_long = p_long * 6894.76;

c_long = linspace(cmin,cmax,npts);

theta_longd = linspace(thetamin,thetamax,npts); % degrees
theta_long = theta_longd*pi/180; %radians

V_long_percent = linspace(Vmin,Vmax,npts);
V_long = V_long_percent*const.V_b;




%Doing calculations again with long matrices for input to plotScatter
[time_thetalong,state_thetalong] = varying_parameters("theta_i",theta_long,int_time);
[time_Plong,state_plong] = varying_parameters("p_0",p_long,int_time);
[time_clong,state_clong] = varying_parameters("c_D",c_long,int_time);
[time_Vlong,state_Vlong] = varying_parameters("V_0w",V_long,int_time);

% Plotting max distance with array of varied parameters for optimization
[max_distanceT,idealTheta] = plotScatter(theta_longd,state_thetalong,'theta_i');
[max_distanceP,idealP] = plotScatter(p_long,state_plong,'p_0');
[max_distanceC,idealC] = plotScatter(c_long,state_clong,'c_D');
[max_distanceV,idealV] = plotScatter(V_long_percent,state_Vlong,'V_0w');

%Reconverting from Pa to PSI
idealPpsi = idealP/6894.76;

%Printing optimal constraints for farthest flight (changing only one
%variable)
fprintf('Best theta for max distance is %2.2f degrees for %2.2f m\n',idealTheta(1),max_distanceT(1));
fprintf('Best c_D for max distance is %2.2f for %2.2f m\n',idealC(1),max_distanceC(1));
fprintf('Best P_0 for max distance is %2.2f psi for %2.2f m\n',idealPpsi(1),max_distanceP(1));
fprintf('Best percent volume of water for max distance is %2.2f%% for %2.2f m\n',100*idealV(1),max_distanceV(1));

disp(' ');

fprintf('Best theta for max height is %2.2f degrees for %2.2f m\n',idealTheta(2),max_distanceT(2));
fprintf('Best c_D for max height is %2.2f for %2.2f m\n',idealC(2),max_distanceC(2));
fprintf('Best P_0 for max height is %2.2f psi for %2.2f m\n',idealPpsi(2),max_distanceP(2));
fprintf('Best percent volume of water for max height is %2.2f%% for %2.2f m\n',100*idealV(2),max_distanceV(2));

disp(' ');

%% Using all of the optimal conditions for maximum distance

%Using the update4Const function to calculate the trajectory with all 4
%best values
constBest = update4Const("theta_i",idealTheta(1)*(pi/180),"c_D",idealC(1),"p_0",idealP(1),"V_0w",idealV(1)*const.V_b);

%Calculating initial conditions for ode45 (some don't change)
vx0 = constBest.v_0*cos(constBest.theta_i);
vz0 = constBest.v_0*sin(constBest.theta_i);

% Creating a column vector of initial condition
initial_conditions = [constBest.x_0 ; vx0 ; constBest.z_0 ; vz0 ; constBest.m_0tot ; constBest.V_0a ; constBest.m_0a];

%A new time span beacause when optimized the rocket takes more than 5
%seconds to fly
int_time2 = [0,6];

%new function handle with the best const struct
f = @(t,y)state_matrix_func(constBest,t,y);

%Calculating the flight using the best conditions for maximum flight
[best_t,best_state] = ode45(f,int_time2,initial_conditions);

%Plotting the best trajectory
figure();
hold on;
plot(best_state(:,1),best_state(:,3),'linewidth',1);
ylim([0,35]);
xlabel('Distance (m)');
ylabel('Height (m)');
grid minor;
title('Best Possible Trajectory with Optimal Conditions');

%% Hitting 85 m

%Changeable Parameters to adjust to test how close to 80m is possible
target_psi = 80;
target_c = 0.4259;
target_psi = target_psi*6894.76;

%Using update 2 const to change 2 parameters in the const struct
constTarget = update2Const("p_0",target_psi,"c_D",target_c);

vx0 = constTarget.v_0*cos(constTarget.theta_i);
vz0 = constTarget.v_0*sin(constTarget.theta_i);

initial_conditions = [constTarget.x_0 ; vx0 ; constTarget.z_0 ; vz0 ; constTarget.m_0tot ; constTarget.V_0a ; constTarget.m_0a];

%Setting up a new function handle 
f_target = @(t,state)state_matrix_func(constTarget,t,state);


%Calculating the flight using the best conditions for maximum flight
[target_t,target_state] = ode45(f_target,int_time,initial_conditions);

% Preallocating
thrust = zeros(length(target_t),1);
stage = zeros(length(target_t),1);

%Calculating the thrust and stage from the ode45 values
for i = 1:length(target_t)
    [~,thrust(i),stage(i)] = state_matrix_func(constTarget,target_t(i),target_state(i,:));
end


%Index where the change occurs 
stage2 = find(stage==2,1);
stage3 = find(stage==3,1);
stage4 = find(stage==4,1);

% The transitions between stages in time 
transition1_time = target_t(stage2);
transition2_time = target_t(stage3);
transition3_time = target_t(stage4);

%The transition between stages in x position
transition1_x = target_state(stage2,1);
transition2_x = target_state(stage3,1);
transition3_x = target_state(stage4,1);

%Plotting the best trajectory
figure();
hold on;
plot(target_state(:,1),target_state(:,3),'linewidth',1);
ylim([0,25]);
xlabel('Distance (m)');
ylabel('Height (m)');
grid minor;
xline(transition1_x);
xline(transition2_x);
xline(transition3_x);
title('Trajectory to hit 85 m');
legend('Trajectory','Stage Transition');


figure();
hold on;
plot(target_t,thrust,'Linewidth',1);
xline(transition1_time);
xline(transition2_time);
xline(transition3_time);
grid minor;
xlabel('Time (s)');
ylabel('Thrust (N)');
title('Thrust vs. Time');
xlim([0,0.5]);
legend('Thrust','Stage Transitions');

%Printing the landing distance to see if the given parameters provided the
%desired trajectory
fprintf('Distance traveled: %2.2f\n',max(target_state(:,1)));
fprintf('Peak thrust is %2.2f N\n',max(thrust));
fprintf('Max height: %2.2f m\n',max(target_state(:,3)));

%% FUNCTIONS 

function const = update4Const(p1,val1,p2,val2,p3,val3,p4,val4)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function is the same as updateConst, just working for 4 different
% parameters. This function works for exactly 4 values. No more, no less.

% INPUTS:

% p(i) is a string, the name of the struct value being changed. In
% order for this function to work, it must be identical to the name of the
% value in struct

% val(i) is the value being changed to. 

% OUTPUTS:

% const is the const struct that will be used in the calculations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    const.g = 9.807; %m/s^2 (gravitational constant)
    const.c_dis = 0.8; %(Discharge constant)
    const.rho_air = 0.961; %kg/m^3 (Density of ambient air)
    const.V_b = 0.002; %m^3 (Volume of bottle)
    const.p_atm = 12.1 * 6894.76; %psia to Pa (pressure of atmosphere)
    const.gamma = 1.4; %unitless (specific heat reatio constant)
    const.rho_w = 1000; %kg/m^3 (density of water)
    const.d_e = 2.1; %cm (diameter of exit)
    const.d_b = 10.5; %cm (diameter of bottle)
    const.R_air = 287; %J/kgK (Gas constant for air)
    const.m_b = 0.15; %kg (mass of bottle)
    const.c_D = 0.48; %Coefficient of drag)
    const.p_0 = 52* 6894.76 + const.p_atm; %psig to Pa (initial pressure in bottle)
    const.V_0w = 0.00095; %m^3 (Initial volume of water)
    const.T_0 = 300; %k (Initial temperature of air)
    const.v_0 = 0.0; %m/s (initial velocity)
    const.theta_i = 42*(pi/180); % degrees to radian (launch angle)
    const.x_0 = 0; %m (initial x position)
    const.z_0 = 0.25; %m (initial z position)
    const.l_s = 0.5; %m (length of launch stand)

    %Reassigning the values
    const.(p1) = val1;
    const.(p2) = val2;
    const.(p3) = val3;
    const.(p4) = val4;

    %Calculating values that may or may not change based on the inputted
    %parameter values
    const.At  =   pi*((const.d_e/2)*0.01)^2; %m^2 (Cross sectional area of throat)
    const.Ab  =   pi*((const.d_b/2)*0.01)^2; %m^2 (Cross sectional area of bottle)

    const.m_0w = const.rho_w * const.V_0w; %kg (initial mass of water)
    rho_0a = const.p_0/(const.R_air*const.T_0); %kg/m^3 (initial density of air)
    const.V_0a = const.V_b - const.V_0w; %m^3 (Initial volume of air)
    const.m_0a = const.V_0a * rho_0a; %kg (initial mass of air)

    const.m_0tot = const.m_b + const.rho_w * (const.V_b - const.V_0a) + (const.p_0*const.V_0a)/(const.R_air*const.T_0);

end
function [state_matrix,F,stage] = state_matrix_func(const,~,state) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function that takes in the current state of the rocket, and has
% equations for calculating the derivatives based on current state

% This function will output the derivatives of each part of the state
% matrix

% state_matrix = d/dt[x;vx;z;vz;mr;V;m_air]

% This function outputs the values of each of these derivatives for ode45,
% and F, the thrust, and stage, an integer value of the stage the rocket is
% in. 

% This function requires input of the constant structure, time, and the
% current state since the derivatives are dependant on the current state. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Constants 
   %These come from the const struct that is inputted

    g = const.g; %m/s^2
    c_dis = const.c_dis;
    rho_air = const.rho_air; %kg/m^3
    V_b = const.V_b; %m^3
    p_atm = const.p_atm; %pa
    gamma = const.gamma;
    rho_w = const.rho_w; %kg/m^3
    R_air = const.R_air; %J/kgK
    c_D = const.c_D;
    p_0 = const.p_0; %pa
    theta_i = const.theta_i; %degree
    z0 = const.z_0; %m
    l_s = const.l_s; %m
    At = const.At; %m^2
    Ab = const.Ab; %m^2
    V_0a = const.V_0a; %m^3
    m_0a = const.m_0a; %kg

    % Lowercase v is velocity, uppercase V is  volume

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Devectorizing the State matrix
    x1 = state(1);
    x2 = state(2);

    z1 = state(3);
    z2 = state(4);

    m_r = state(5);
    V_air = state(6);
    m_air = state(7);

     %Calculating distance and theta, to see if the bottle is of the launch
     %stand, and locking theta if it is. 
    distance = sqrt(x1^2+(z1-z0)^2);

    if distance < l_s
        theta = theta_i;
    else
        theta = atan(z2/x2);
    end
    
    %Calculating dynamic pressure and drag at any point in the integration
    v_bottle = sqrt(x2^2+z2^2);
    q = 0.5 * rho_air * (v_bottle^2);
    drag = q*Ab*c_D;

    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Stages of flight


    if V_air < V_b  % Stage 1: Water expulsion

        %Calculating pressure of air
        p = p_0*(V_0a/V_air)^gamma;
        

        %Calculating state of exhaust
        v_e = sqrt(2*(p-p_atm)/rho_w); %Velocity of the exhausted water

        %calulating mass flow rate of water
        mdot_w = c_dis * rho_w * At * v_e;

        %Calculating thrust
        F = mdot_w * v_e;  %Thrust created by exhausted water
        
        %Calculating the time rate of change of water
        dVdt = c_dis * At * v_e;

        %Finding net forces in x and z directions
        FnetZ = F*sin(theta)-drag*sin(theta)-m_r*g;
        FnetX = F*cos(theta)-drag*cos(theta);

        ddt1 = x2; %x Velocity
        ddt2 = FnetX/m_r; %Thrust - drag
        ddt3 = z2; %Z velocity
        ddt4 = FnetZ/m_r;  %Thrust-drag-gravity
        ddt5 = -mdot_w; % Current total mass change in kg/s
        ddt6 = dVdt; %Change in volume of air
        ddt7 = 0; %Change in mass of air
        %disp('Water expulsion');
        stage = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    elseif V_air >= V_b   % Stages 2-4 


        %Calculating P end and T end for calculating pressure later
        p_end = p_0*(V_0a/V_b)^gamma;
        

        %Calculating density, Temperature, and pressure of air
        rho = m_air/V_b;
        p = p_end*(m_air/m_0a)^gamma;
        T = p/(rho*R_air);
        
        %these are the actual current temperature and density     

        %Determining if the flow is choked to find v_e
        p_star = p*(2/(gamma+1))^(gamma/(gamma-1));

        if p_star > p_atm
            T_e = (2/(gamma+1))*T;
            rho_e = p_star/(R_air*T_e);
            v_e = sqrt(gamma*R_air*T_e);
            p_e = p_star;
            
        else
            M = sqrt(((p/p_atm)^((gamma-1)/gamma)-1)*(2/(gamma-1)));
            T_e = T*(1+(gamma-1)/2*M^2)^-1;
            rho_e = p_atm/(R_air*T_e);
            p_e = p_atm;
            
            v_e = M*sqrt(gamma*R_air*T_e);
            
            
        end


        %Calculating thrust
        m_dota = c_dis*rho_e*At*v_e;
        F = m_dota*v_e + (p_e-p_atm)*At; 

        %Calculating change in total mass
        mdot_r = -c_dis*rho_e*At*v_e;

        FnetZ = F*sin(theta)-drag*sin(theta)-m_r*g;
        FnetX = F*cos(theta)-drag*cos(theta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if p>p_atm %Stage 2: Air Expulsion

            ddt1 = x2; %x Velocity
            ddt2 = FnetX/m_r; %Thrust - drag
            ddt3 = z2; %Z velocity
            ddt4 = FnetZ/m_r;  %Thrust-drag-gravity
            ddt5 = mdot_r;
            ddt6 = 0;
            ddt7 = mdot_r;
            %disp('Air Expulsion');
            stage = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif p <= p_atm && z1 > 0  %Stage 3: Ballistic Stage

            ddt1 = x2;
            ddt2 = -cos(theta)*drag/m_r;
            ddt3 = z2;
            ddt4 = -sin(theta)*drag/m_r-g;
            ddt5 = 0;
            ddt6 = 0;
            ddt7 = 0;
            %disp('Ballistic');
            stage = 3;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif z1<=0   %Stage 4: Landed (no change in any variable)

             ddt1 = 0;
             ddt2 = 0;
             ddt3 = 0;
             ddt4 = 0;
             ddt5 = 0;
             ddt6 = 0;
             ddt7 = 0;
             %disp('Landing');
             stage = 4;

        end

    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    else  %Error phase just in case 
            ddt1 = 0;
            ddt2 = 0;
            ddt3 = 0;
            ddt4 = 0;
            ddt5 = 0;
            ddt6 = 0;
            ddt7 = 0;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Creating the state matrix with the parts to integrate 
%This will be outputted by the function and will be put into ode45

    state_matrix = [ddt1;ddt2;ddt3;ddt4;ddt5;ddt6;ddt7];

end
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
function const = getConst()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function creates a structure of initial conditions and relevant
% constants for the main function. This truncates the state_matrix_func
% inputs by taking in only this structure. 

% No inputs are required for the function. 
% Output is a structure with multiple different entries

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    const.g = 9.807; %m/s^2 (gravitational constant)
    const.c_dis = 0.8; %(Discharge constant)
    const.rho_air = 0.961; %kg/m^3 (Density of ambient air)
    const.V_b = 0.002; %m^3 (Volume of bottle)
    const.p_atm = 12.1 * 6894.76; %psia to Pa (pressure of atmosphere)
    const.gamma = 1.4; %unitless (specific heat reatio constant)
    const.rho_w = 1000; %kg/m^3 (density of water)
    const.d_e = 2.1; %cm (diameter of exit)
    const.d_b = 10.5; %cm (diameter of bottle)
    const.R_air = 287; %J/kgK (Gas constant for air)
    const.m_b = 0.15; %kg (mass of bottle)
    const.c_D = 0.48; %Coefficient of drag)
    const.p_0 = 52* 6894.76 + const.p_atm; %psig to Pa (initial pressure in bottle) 52* 6894.76 + const.p_atm
    const.V_0w = 0.00095; %m^3 (Initial volume of water)
    const.T_0 = 300; %k (Initial temperature of air)
    const.v_0 = 0.0; %m/s (initial velocity)
    const.theta_i = 42*(pi/180); % degrees to radian (launch angle)
    const.x_0 = 0; %m (initial x position)
    const.z_0 = 0.25; %m (initial z position)
    const.l_s = 0.5; %m (length of launch stand)
   
    
    % Calculating other necessary constants and initial conditions
    const.At  =   pi*((const.d_e/2)*0.01)^2; %m^2 (Cross sectional area of throat)
    const.Ab  =   pi*((const.d_b/2)*0.01)^2; %m^2 (Cross sectional area of bottle)

    const.m_0w = const.rho_w * const.V_0w; %kg (initial mass of water)
    rho_0a = const.p_0/(const.R_air*const.T_0); %kg/m^3 (initial density of air)
    const.V_0a = const.V_b - const.V_0w; %m^3 (Initial volume of air)
    const.m_0a = const.V_0a * rho_0a; %kg (initial mass of air)

    const.m_0tot = const.m_b + const.rho_w * (const.V_b - const.V_0a) + (const.p_0*const.V_0a)/(const.R_air*const.T_0);
end
function plotFun(state,parameter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function takes in cells of time and state from ode45 and plots all
% of the struct rows on similar graphs

% No return is necessary, plots will just be created


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot of trajectory 
    figure()
    hold on;
    
    %For loop to plot all of the 8 different lines
    for i = 1:8
        plot(state{i}(:,1),state{i}(:,3),'Linewidth',1);
    end

    grid on;
    xlabel ('Distance (m) ');
    ylabel ( 'Distance (m)');


    text = parameter;

    title(sprintf('Trajectory based on changing %s',text));
    legend(sprintf('%s 1',text),sprintf('%s 2',text), ...
        sprintf('%s 3',text),sprintf('%s 4',text), ...
        sprintf('%s 5',text),sprintf('%s 6',text), ...
        sprintf('%s 7',text),sprintf('%s 8',text));

    hold off;

end
function const = update2Const(p1,val1,p2,val2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function is the same as updateConst, just working for 4 different
% parameters. This function works for exactly 2 values. No more, no less.

% INPUTS:

% p(i) is a string, the name of the struct value being changed. In
% order for this function to work, it must be identical to the name of the
% value in struct

% val(i) is the value being changed to. 

% OUTPUTS:

% const is the const struct that will be used in the calculations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    const.g = 9.807; %m/s^2 (gravitational constant)
    const.c_dis = 0.8; %(Discharge constant)
    const.rho_air = 0.961; %kg/m^3 (Density of ambient air)
    const.V_b = 0.002; %m^3 (Volume of bottle)
    const.p_atm = 12.1 * 6894.76; %psia to Pa (pressure of atmosphere)
    const.gamma = 1.4; %unitless (specific heat reatio constant)
    const.rho_w = 1000; %kg/m^3 (density of water)
    const.d_e = 2.1; %cm (diameter of exit)
    const.d_b = 10.5; %cm (diameter of bottle)
    const.R_air = 287; %J/kgK (Gas constant for air)
    const.m_b = 0.15; %kg (mass of bottle)
    const.c_D = 0.48; %Coefficient of drag)
    const.p_0 = 52* 6894.76 + const.p_atm; %psig to Pa (initial pressure in bottle)
    const.V_0w = 0.00095; %m^3 (Initial volume of water)
    const.T_0 = 300; %k (Initial temperature of air)
    const.v_0 = 0.0; %m/s (initial velocity)
    const.theta_i = 42*(pi/180); % degrees to radian (launch angle)
    const.x_0 = 0; %m (initial x position)
    const.z_0 = 0.25; %m (initial z position)
    const.l_s = 0.5; %m (length of launch stand)

    %Reassigning the values
    const.(p1) = val1;
    const.(p2) = val2;


    %Calculating values that may or may not change based on the inputted
    %parameter values
    const.At  =   pi*((const.d_e/2)*0.01)^2; %m^2 (Cross sectional area of throat)
    const.Ab  =   pi*((const.d_b/2)*0.01)^2; %m^2 (Cross sectional area of bottle)

    const.m_0w = const.rho_w * const.V_0w; %kg (initial mass of water)
    rho_0a = const.p_0/(const.R_air*const.T_0); %kg/m^3 (initial density of air)
    const.V_0a = const.V_b - const.V_0w; %m^3 (Initial volume of air)
    const.m_0a = const.V_0a * rho_0a; %kg (initial mass of air)

    const.m_0tot = const.m_b + const.rho_w * (const.V_b - const.V_0a) + (const.p_0*const.V_0a)/(const.R_air*const.T_0);

end
function [tout,states] = varying_parameters(parameter,test_parameter,tspan)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is made to run multiple calculations of different initial
% parameters.
% A for loop will call ode45 with the const struct changing for each new
% parameter value being changed. Only one parameter can be changed at a
% time, to change multiple different ones, this function must be called
% again. 

% INPUTS:
% Parameter is the desired parameter being changed, as a string. This
% should be the name in the const struct (from getConst()) that is being
% changed. 

% const is the struct that is used to calculate the trajectory with given
% initial parameters. 

% test_parameter is a matrix of values that the changing parameter is
% spanning 

% tspan is the required integration time span for ode45

% OUTPUTS:
% tout is a cell of length(test_parameter), where each entry is a matrix of
% time values calculated from ode45.

% states is a cell of the same length where each entry is a matrix of the
% state matrix calculated by ode45.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Preallocation
tout = cell(length(test_parameter),1);
states = cell(length(test_parameter),1);


    for i = 1:length(test_parameter)
        
        % Changing the parameter inside the const struct to the new value
        const = updateConst(parameter,test_parameter(i));
        
        
        %Calculating initial conditions for ode45 (some don't change)
        vx0 = const.v_0*cos(const.theta_i);
        vz0 = const.v_0*sin(const.theta_i);

        % Creating a column vector of initial conditions
        initial_conditions = [const.x_0 ; vx0 ; const.z_0 ; vz0 ; const.m_0tot ; const.V_0a ; const.m_0a];

        %Function handle, passing in the new const struct
        f = @(t,state)state_matrix_func(const,t,state);

        %calculating the output values with the given function 
        [tout{i},states{i}] = ode45(f,tspan,initial_conditions);
        
    end
    

end
function const = updateConst(parameter, test_parameter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function is designed to change parameters of the const struct based 
% inputted new values

% INPUTS:

% paramater is a string, the name of the struct value being changed. In
% order for this function to work, it must be identical to the name of the
% value in struct

% test_parameter is the value being changed to. 

% OUTPUTS:

% const is the const struct that will be used in the calculations


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    const.g = 9.807; %m/s^2 (gravitational constant)
    const.c_dis = 0.8; %(Discharge constant)
    const.rho_air = 0.961; %kg/m^3 (Density of ambient air)
    const.V_b = 0.002; %m^3 (Volume of bottle)
    const.p_atm = 12.1 * 6894.76; %psia to Pa (pressure of atmosphere)
    const.gamma = 1.4; %unitless (specific heat reatio constant)
    const.rho_w = 1000; %kg/m^3 (density of water)
    const.d_e = 2.1; %cm (diameter of exit)
    const.d_b = 10.5; %cm (diameter of bottle)
    const.R_air = 287; %J/kgK (Gas constant for air)
    const.m_b = 0.15; %kg (mass of bottle)
    const.c_D = 0.48; %Coefficient of drag)
    const.p_0 = 52* 6894.76 + const.p_atm; %psig to Pa (initial pressure in bottle) 
    const.V_0w = 0.00095; %m^3 (Initial volume of water)
    const.T_0 = 300; %k (Initial temperature of air)
    const.v_0 = 0.0; %m/s (initial velocity)
    const.theta_i = 42*(pi/180); % degrees to radian (launch angle)
    const.x_0 = 0; %m (initial x position)
    const.z_0 = 0.25; %m (initial z position)
    const.l_s = 0.5; %m (length of launch stand)

    % Changing the desired parameter value
    const.(parameter) = test_parameter;

    %Finding all the values that have to be calculated, after the new value
    %has been changed out
    
    const.At  =   pi*((const.d_e/2)*0.01)^2; %m^2 (Cross sectional area of throat)
    const.Ab  =   pi*((const.d_b/2)*0.01)^2; %m^2 (Cross sectional area of bottle)

    const.m_0w = const.rho_w * const.V_0w; %kg (initial mass of water)
    rho_0a = const.p_0/(const.R_air*const.T_0); %kg/m^3 (initial density of air)
    const.V_0a = const.V_b - const.V_0w; %m^3 (Initial volume of air)
    const.m_0a = const.V_0a * rho_0a; %kg (initial mass of air)

    const.m_0tot = const.m_b + const.rho_w * (const.V_b - const.V_0a) + (const.p_0*const.V_0a)/(const.R_air*const.T_0);

end

