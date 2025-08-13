function R321 = RotationMatrix321(attitude)
% Standard 321 Euler Angle Rotation
% If xI is inertial coordinates of vector
% then XB = R321 * XI are body coordinates
%
phi = attitude(1);
theta = attitude(2);
psi = attitude(3);

R321 = [cos(theta)*cos(psi), cos(theta)*sin(psi), -sin(theta); 
    sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi), sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi), sin(phi)*cos(theta); 
    cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi), cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi), cos(phi)*cos(theta)];
