function [e,c_L,c_Di,delta] = PLLT2(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)

% INPUTS: 
% b - span
% a0_t - lift slope at the wing tip
% a0_r - lift slope at the wing root
% c_t - chord length at the tip
% c_r - chord length at root
% aero_t - induced AoA at wing tip
% aero_r - induced AoA at wing root
% geo_t - geometric AoA at wing tip
% geo_r - geometric AoA at wing root
% N - number of Fourier Coefficients (terms)


% OUTPUTS:
% e - span efficiency factor
% c_L - coefficient of lift
% c_Di - induced drag coefficient


%% Calculating the geometric parameters as functions to use later
S = 2*(0.5*(c_t + c_r)*b/2);
AR = b^2/S;
c = @(y) (c_t-c_r)/(-b/2)*y + c_r; % Assuming linear change
a0 = @(y) (a0_t-a0_r)/(-b/2)*y + a0_r;
aero = @(y) (aero_t-aero_r)/(-b/2)*y + aero_r;
geo = @(y) (geo_t-geo_r)/(-b/2)*y + geo_r;

% Control points
theta0 = zeros(N,1);
for i = 1:N
    theta0(i) = i*pi/(2*N);
end

y0 = -b/2*cos(theta0);


%% Solving system of equations to find the Fourier Coefficients
M = zeros(N,N); %Matrix for linear system of equations

% Assigning the entries of B
for i = 1:N
    for j = 1:N
        M(i,j) = 4*b/(a0(y0(i))*c(y0(i)))*sin((2*j-1)*theta0(i)) + (2*j-1)*sin((2*j-1)*theta0(i))/sin(theta0(i));
    end
end

% Assigning the entries of b
rhs = zeros(N,1);
for i = 1:N
    rhs(i) = -aero(y0(i)) + geo(y0(i));
end

% Solving for the Fourier Coefficients
%A = M\rhs;
A = inv(M)*rhs;
% A1 should be 0.01261

%% Solving for the desired quantities
c_L = A(1)*pi*AR;
delta = 0;

for i = 2:N
    delta = (2*i-1)*(A(i)/A(1))^2 + delta;
end
e = 1/(1+delta);
c_Di = c_L^2/(pi*e*AR);

end