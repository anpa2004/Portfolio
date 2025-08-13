function xdot = AircraftEOMDoublet(time, aircraft_state, aircraft_surfaces,doublet_size, doublet_time, wind_inertial, aircraft_parameters)
% This function is used instead of Aircraft EOM, where it uses AircraftEOM
% but changes the elevator deflection depending on the time. This is used
% to simulate Phugoid mode


de0 = aircraft_surfaces(1);
% Changing the deflection based on the time and size inputted
if  time<doublet_time
    de = de0+doublet_size;
    
end
if (time>doublet_time) && (time<2*doublet_time)
    de = de0-doublet_size;
    
end
if (time>2*doublet_time)
    de = de0;
    
end
aircraft_surfaces(1) = de;

% Calculating the state derivatives given the deflection
xdot = AircraftEOM(time, aircraft_state, aircraft_surfaces,wind_inertial, aircraft_parameters);

