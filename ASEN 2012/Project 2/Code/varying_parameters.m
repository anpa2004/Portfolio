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