function [wind_angles] = WindAnglesFromVelocityBody(air_rel_vel_body)
%This function takes the inputs from the relative air velocity, and uses
%them to caluclate the magnitude of the airspeed and the respective angles
%alpha and beta.

u = air_rel_vel_body(1);
v = air_rel_vel_body(2);
w = air_rel_vel_body(3);

V = sqrt(u^2 + v^2 + w^2);
 beta = asin(v/V);
% alpha = acos(u/(V*cos(beta)));
alpha = atan2(w,u);

wind_angles = [V;beta;alpha];

end

