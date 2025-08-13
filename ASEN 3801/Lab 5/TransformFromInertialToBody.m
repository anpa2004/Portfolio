function vector_body = TransformFromInertialToBody(vector_inertial, euler_angles)
vector_body = RotationMatrix321(euler_angles)*vector_inertial;


% function wind_body = TransformFromInertialToBody(wind_inertial, aircraft_state)
% %The inputs are the aircraft state vector inputs are the entries 4:6, aka
% %the Euler angles, and the wind vector in the inertial frame. 
% 
 
