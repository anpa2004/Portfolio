function [ddt] = EOM(t,state,const)
    
    %Taking const values out of struct for ease
    g = const.g;
    c_dis = const.c_dis;
    rhoair = const.rhoair;
    V_b = const.V_b;
    p_atm = const.patm;
    gamma = const.gamma;
    rhow = const.rhow;
    d_e = const.d_e;
    d_b = const.d_b;
    R = const.R;
    c_D = const.c_D;
    p_0 = const.p_0;
    theta0 = const.theta;
    l_s= const.l_s;
    r0 = const.r0;
    A_b = pi*(d_b/2)^2;
    A_t = pi*(d_e/2)^2;
    V_0a = const.V_0a;
    m_a=const.m_a;

    
    %Devectorizing the state vector for ease later
    r = [state(1);state(3)];
    v = [state(2);state(4)];
    m_r = state(5);
    V_air = state(6);
    m_air = state(7);

    %Dealing with the launch stand
    x_s = l_s*cos(theta0);
    z_s = l_s*sin(theta0);
    
    %This is the heading the rocket will have if it is still on the stand
    h_stand = [x_s;z_s];
    h_stand = h_stand/norm(h_stand);

    %How far the rocket has travelled from the launch point
    displacement = r-r0;
    
    %If rocket is still on the launch stand
    if norm(displacement)<=l_s
        v = norm(v)*h_stand; %Speed*direction
    end
    
    heading = v/norm(v);
    %Calculating pressure drag, drag vector
    q = 0.5*rhoair*v.^2;
    
    drag = q*c_D*A_b;

    F_w = [0;m_r*g];


    %% Phase 1: Water Expulsion
    if v_air<V_b
        %Pressure of air at any given time
        p = p_0*(V_0a/V_air)^gamma;
        
        %Exhaust velocity
        v_e = sqrt(2*(p-p_atm)/rhow);
        
        %calculating thrust
        F = 2*c_dis*A_t*(p-p_atm); %scalar of force
        thrust = F*heading; %vector
        
        %Time rate of change of volume of air
        dVdt = c_dis*A_t*v_e;
        
        %Mass flow rate of water
        m_dotW = c_dis*rhow*A_t*v_e;
        m_dotR = - m_dotW;


        %Calculating Forces
        

        %Calculating Net Force, acceleration
        F_net = F_w+thrust-drag;

        a = F_net/m_r;

        %Assigning the derivative values to the values of the EOM return
        ddt1 = v(1);
        ddt2 = a(1);
        ddt3 = v(2);
        ddt4 = a(2);
        ddt5 = m_dotR;
        ddt6 = dVdt;
        ddt7 = 0;


    else
        p_end = p_0*(V_0a/V_b)^gamma;

        p = p_end*(m_air/m_a);

        % Phase 2: Air Exhaustion
        if p>p_atm

            %Determining if flow is choked                       
            p_star = p*(2/(gamma+1))^(gamma/(gamma-1));
            rhoA = m_air/V_b;
            T = p/(rhoA*R);

            if p_star>p_atm
                M_e = 1;
                T_e = (2/(1+gamma))*T;
                p_e = p_star;
                rho_e = p_e/(R*T_e);
              
            else
                M_e = sqrt((2*((p/p_atm)^((gamma-1)/gamma)-1))/(gamma-1));
                T_e = T*(1+((gamma-1)/2)*M_e^2)^-1;
                p_e = p_atm;
                rho_e = p_e/(R*T_e);
            end

            %Calculate velocity of exhaust
            v_e = M_e*sqrt(gamma*R*T_e);
            
            % Mass flow rate of air
            m_dotA = c_dis*rho_e*A_t*v_e;

            % Force from thrust
            F = m_dotA*v_e + (p_e-p_atm)*A_t;
            
            % Rate of chane of total mass
            m_dotR = - m_dotA;
            
            %Vector of thrust, acceleration.
            thrust = F*heading;
            F_net = F_w+thrust-drag;

            a = F_net/m_r;

            ddt1 = v(1);
            ddt2 = a(1);
            ddt3 = v(2);
            ddt4 = a(2);
            ddt5 = -m_dotA;
            ddt6 = 0;
            ddt7 = m_dotR;
        else
        %Ballistic Phase
            %Calculating net force (drag and weight only)
            F_net = F_w - drag;
            a = F_net/m_r;

            ddt1 = v(1);
            ddt2 = a(1);
            ddt3 = v(2);
            ddt4 = a(2);
            ddt5 = 0;
            ddt6 = 0;
            ddt7 = 0;
            
        end
    end
    %% Creating Derivative vector
    ddt= [ddt1;ddt2;ddt3;ddt4;ddt5;ddt6;ddt7];

end