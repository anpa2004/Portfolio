function var_dot = QuadrotorEOM(t,var,g,m,I,d,km,nu,mu,motor_forces)

phi = var(4);
theta = var(5);
psi = var(6);

uE = var(7);
vE = var(8);
wE = var(9);

p = var(10);
q = var(11);
r = var(12);

Ix = I(1,1);
Iy = I(2,2);
Iz = I(3,3);

Va = sqrt(uE^2+vE^2+wE^2);

f1= motor_forces(1);
f2 = motor_forces(2);
f3 = motor_forces(3);
f4 = motor_forces(4);

L = -mu*sqrt(p^2+q^2+r^2)*p;
M = -mu*sqrt(p^2+q^2+r^2)*q;
N = -mu*sqrt(p^2+q^2+r^2)*r;

A = [-1 -1 -1 -1; -d/sqrt(2) -d/sqrt(2) d/sqrt(2) d/sqrt(2); d/sqrt(2) -d/sqrt(2) -d/sqrt(2) d/sqrt(2); km -km km -km];
v0 = [f1; f2; f3; f4];

v = A*v0;

% Perturbation to determine natural stability

% if t>2&&t<2.2
%     
%    theta = 0.5; 
% 
% end

Zc = v(1);
Lc = v(2);
Mc = v(3);
Nc = v(4);



% Vector EOM 1
A1 = [cos(theta)*cos(psi) sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi);
    cos(theta)*sin(psi) sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi);
    -sin(theta) sin(phi)*cos(theta) cos(phi)*cos(theta)];
v1 = [uE;vE;wE];
v2 = A1*v1;
var_dot(1)=v2(1);
var_dot(2) = v2(2);
var_dot(3) = v2(3);

% Vector EOM 2
A2 = [1 sin(phi)*tan(theta) cos(phi)*tan(theta);0 cos(phi) -sin(phi); 0 sin(phi)*sec(theta) cos(phi)*sec(theta)];
v3 = [p;q;r];
v4 = A2*v3;



var_dot(4)=v4(1);
var_dot(5) = v4(2);
var_dot(6) = v4(3);

% Vector EOM 3
var_dot(7) = r*vE-q*wE -g*sin(theta) -nu*Va*uE/m;
var_dot(8) = p*wE - r*uE + g*cos(theta)*sin(phi) -nu*Va/m*vE;
var_dot(9) = q*uE - p*vE +g*cos(theta)*cos(phi)-nu*Va/m*wE + Zc/m;

% Vector EOM 4
var_dot(10) = (Iy-Iz)/Ix*q*r + L/Ix + Lc/Ix;
var_dot(11) = (Iz-Ix)/Iy*p*r + M/Iy + Mc/Iy;
var_dot(12) = (Ix-Iy)/(Iz)*p*q + N/Iz+Nc/Iz;


var_dot = var_dot';