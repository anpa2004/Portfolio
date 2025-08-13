function ddt = ddt_fun(t,state,const)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function calculates the derivatives of the balloon's state for
% ode45 to integrate

%INPUTS: 
% t is the time value (no functions are dependant on t so this is unused)
% State is the current state, a column vector of x,vx,z,vz.
% const is a struct of all initial relevant constants. 

%OUTPUTS:

% dzdt is a column vector of the derivatives of state (d/dt[x;vx;z;vz])
% that will be integrated to model the flight of the balloon

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Devectorizing the y vector for simplicity later
    vx = state(2);
    z = state(3);
    vz = state(4);
  
    %Taking constants and scalars from the const struct for simplicity
    %later
    mtot = const.mass_tot; %kg
    mPay = const.mass_payload; %kg
    rhoA = const.rhoAir; %kg/m^3
    g = const.g; %m/s^2
    C1 = const.cd1;
    C2 = const.cd2;
    C3 = const.cd3;
    A = const.A; %m^2
    V = const.V; %m^3 
    pop = const.hpop; %m (altitude of pop)
    chute = const.hchute; %m (altitude of chute deployment)


    %velocities of wind and balloon
    v_wind = [const.v_wind;0];
    v_balloon = [vx;vz];

    %relative velocity of balloon
    v_rel = v_balloon - v_wind;

    %Heading 
    h = v_rel/norm(v_rel);

    %If statements to determine which portion of the flight the balloon is
    %in, based on altitude and velocity direction
    if (z<pop && vz>0)
       %Balloon is rising under bouyant force
       D = .5*C1*A*rhoA*norm(v_rel.^2);
       F_drag = -h*D;
       F_buoy = [0;rhoA*V*g];
       F_weight = [0;-mtot*g];
       acceleration = (1/mtot)*(F_drag+F_weight+F_buoy);
       %disp('stage 1');
    
    elseif (z>pop)
       %Balloon has popped, is still rising because it had momentum
       D = .5*C2*A*rhoA*norm(v_rel.^2);
       F_drag = -h*D;
       F_buoy = [0;0];
       F_weight = [0;-mPay*g];

       acceleration = (1/mPay)*(F_drag+F_weight+F_buoy);
       %disp('stage 2');

    elseif(z<pop && z>chute && vz<0)
        %Payload is falling but the chute has not opened
       D = .5*C2*A*rhoA*norm(v_rel.^2);
       F_drag = -h*D;
       F_buoy = [0;0];
       F_weight = [0;-mPay*g];
       acceleration = (1/mPay)*(F_drag+F_weight+F_buoy); 
       %disp('stage 3');

    else
        %Payload is falling under the chute
       D = .5*C3*A*rhoA*norm(v_rel.^2);
       F_drag = -h*D;
       F_buoy = [0;0];
       F_weight = [0;-mPay*g];
       acceleration = (1/mPay)*(F_drag+F_weight+F_buoy);
       %disp('stage 4');
    end


    %Creating a column vector of the derivatives of overrightarrow{y}
    ddt = [vx;acceleration(1);vz;acceleration(2)];


end