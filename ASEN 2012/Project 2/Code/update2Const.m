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