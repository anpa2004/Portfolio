function [idealConst] = target_function(p1,p1lim,p2,p2lim,p3,p3lim,p4,p4lim,const,target,tolerance)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function is designed to automatically adjust the initial conditions
% in order to hit a target difference with a bottle rocket.

% INPUTS:

% p1-p4 are strings that are being changed (p_0,V_0w,theta_i,c_D)
% const is the constant struct that is allowed to change
% Taget is the desired distance in meters (probably should be
% 50<target<120)

% tolerance is the desired proximity to the target (should not be 0, should
% be greater than 0)

% OUTPUTS:

% idealConst is the const struct with the ideal values updated inside 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if tolerance == 0
        disp('Error, tolerance cannot equal 0');
    end

    vx0 = const.v_0*cos(const.theta_i);
    vz0 = const.v_0*sin(const.theta_i);
    
    % Creating a column vector of initial conditions
    initial_conditions = [const.x_0 ; vx0 ; const.z_0 ; vz0 ; const.m_0tot ; const.V_0a ; const.m_0a];

    int_time = [0 5];

    %Function handle of state matrix function
    int_fun = @(t,state) state_matrix_func(const,t,state);

    %Integrating using ode45
    [~,y] = ode45(int_fun,int_time,initial_conditions);

    diff1 = target-max(y(:,1));

    percent = 0.05;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while abs(diff1)>tolerance
        %Increasing p1 by 5%
        if const.(p1)<p1lim(2) && const.(p1)>p1lim(1)
            newVal1 = (1+percent)*const.(p1);
            const = updateConst(p1,newVal1);
            initial_conditions = [const.x_0 ; vx0 ; const.z_0 ; vz0 ; const.m_0tot ; const.V_0a ; const.m_0a];
            [~,y] = ode45(int_fun,int_time,initial_conditions);

            diff2 = abs(target - max(y(:,1)));

            if diff1<diff2
                %Change is bad
                newVal1 = (1-percent)*cpnst.(p1);
                const = updateConst(p1,newVal1);
                initial_conditions = [const.x_0 ; vx0 ; const.z_0 ; vz0 ; const.m_0tot ; const.V_0a ; const.m_0a];
                [~,y] = ode45(int_fun,int_time,initial_conditions);

                diff3 = abs(target-max(y(:,1)));

                %If neither difference is closer to tolerance, it might be
                %near a max, decrease the change parameter
                if diff3>diff1
                    percent = percent-0.01;

                    if percent<0
                        break
                    end

                end
            end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if const.(p2)<p2lim(2) && const.(p2)>p2lim(1)
            newVal2 = (1+percent)*const.(p2);
            const = updateConst(p2,newVal2);
            initial_conditions = [const.x_0 ; vx0 ; const.z_0 ; vz0 ; const.m_0tot ; const.V_0a ; const.m_0a];
            [~,y] = ode45(int_fun,int_time,initial_conditions);

            diff2 = abs(target - max(y(:,1)));

            if diff1<diff2
                %Change is bad
                newVal2 = (1-percent)*cpnst.(p2);
                const = updateConst(p2,newVal2);
                initial_conditions = [const.x_0 ; vx0 ; const.z_0 ; vz0 ; const.m_0tot ; const.V_0a ; const.m_0a];
                [~,y] = ode45(int_fun,int_time,initial_conditions);

                diff3 = abs(target-max(y(:,1)));

                %If neither difference is closer to tolerance, it might be
                %near a max, decrease the change parameter
                if diff3>diff1
                    percent = percent-0.01;

                    if percent<0
                        break
                    end

                end
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if const.(p3)<p3lim(2) && const.(p3)>p3lim(1)
            newVal3 = (1+percent)*const.(p3);
            const = updateConst(p1,newVal3);
            initial_conditions = [const.x_0 ; vx0 ; const.z_0 ; vz0 ; const.m_0tot ; const.V_0a ; const.m_0a];
            [~,y] = ode45(int_fun,int_time,initial_conditions);

            diff2 = abs(target - max(y(:,1)));

            if diff1<diff2
                %Change is bad
                newVal3 = (1-percent)*const.(p3);
                const = updateConst(p3,newVal3);
                initial_conditions = [const.x_0 ; vx0 ; const.z_0 ; vz0 ; const.m_0tot ; const.V_0a ; const.m_0a];
                [~,y] = ode45(int_fun,int_time,initial_conditions);

                diff3 = abs(target-max(y(:,1)));

                %If neither difference is closer to tolerance, it might be
                %near a max, decrease the change parameter
                if diff3>diff1
                    percent = percent-0.01;

                    if percent<0
                        break
                    end

                end
            end
        end
        end



    end
        

    



end