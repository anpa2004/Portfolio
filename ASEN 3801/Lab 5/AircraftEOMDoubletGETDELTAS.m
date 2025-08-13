function aircraft_surfaces = AircraftEOMDoubletGETDELTAS(time, aircraft_state, aircraft_surfaces,doublet_size, doublet_time, wind_inertial, aircraft_parameters)
% This function is iterated through to get the values of the doublets to
% plot since ode45 doesn't return the value of the doublets.

de0 = aircraft_surfaces(1);
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


