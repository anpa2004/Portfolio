function Cd0 = profile_drag_coeff(ct,cr,cdt,cdr,b,S,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the profile drag coefficient for an airfoil
% given the profile drag and the root and tip, assuming trapezoidal wings.
% It also assumes that cd varies linearly across the wingspan between cdt
% and cdr

% INPUTS:
% ct - chord length of the tip
% cr - chord length at the root
% cdt - coefficent of profile drag at the tip
% crt - coefficient of profile drag at the root
% b - wingspan
% S - planform area
% n - number of integration points

% OUTPUTS:
% Cd0 - profile drag coefficient

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Function handles for c(y) and cd(y)
c = @(y) (ct-cr)/(-b/2)*abs(y) + cr;
cd =  @(y) (cdt-cdr)/(-b/2)*abs(y) + cdr;

% Points to evaluate at 
y = linspace(-b/2,b/2,n);

% Integrand for trapz to work with
integrand = cd(y).*c(y);

% Profile drag coefficient
% Cd0 = trapz(integrand,2)/S;
Cd0 = trapz(y,integrand)/S;



