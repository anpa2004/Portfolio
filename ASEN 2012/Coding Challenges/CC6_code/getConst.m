function const = getConst()
    % This is a void function which can be used in the main code to set a
    % variable to the struct this creates
    const.mass_payload = 450; %kg
    const.r = 17; %m
    const.rhoAir = 1.225; %kg/m^3
    const.rhoHe = 0.1786; %kg/m^3
    const.g = 9.81; %m/s^2
    const.V = (4/3)*pi*(const.r^3); %m^3
    const.A = pi*const.r^2; %m^2
    const.mass_gas = const.V*const.rhoHe; %kg
    const.mass_tot = const.mass_gas+const.mass_payload; %kg
    const.cd1 = 0.5;
    const.cd2 = 0;
    const.cd3 = 1.03;
    const.hpop = 2100; %m
    const.hchute = 1950; %m
    const.hBoulder = 1624; %m
    const.x0 = 0; %m
    const.init_velocityZ = 2; %m/s
    const.init_velocityX = 0; %m/s
    const.v_wind = 8; %m/s
    
end