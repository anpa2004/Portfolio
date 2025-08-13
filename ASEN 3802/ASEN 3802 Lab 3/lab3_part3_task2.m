close all; clear; clc;

%%% PART 3 TASK 2 %%%%%%%%%%%%%%%%%%
cr = 5 + 4/12; %ft
ct = 7 + 7/12; %ft
b = 35 + 10/12; %ft
S = 2*(0.5*(ct + cr)*b/2);

alpha = linspace(-6,10,100);


b = 35 + 10/12; %ft
c_r = 5 + 4/12; %ft
c_t = 7 + 7/12; %ft
a0_t = 0.1185 * 180/pi;
a0_r = 0.1186 * 180/pi;
aero_t = 0;
aero_r = 0;
nums_root = [2,4,12];
nums_tip = [0,0,12];
geo_root = 0;
geo_tip = geo_root + 2*pi/180; %rad
N = 20;
n = 10;

% Fitting to be able to model the cd by alpha
cl = [-0.8,0,.85];
cd = [0.018,0.01,0.0195];

% vector of C_L values to use, and associated fit output
cL_eval_tip = linspace(-0.8,0.85,50);
[P,~] = polyfit(cl,cd,2);
cd_eval_tip = polyval(P,cL_eval_tip);

clr = [-0.8,0.2,1.1];
cdr = [0.0218,.009,0.0225];

% vector of C_L values to use, and associated fit output
cL_eval_root = linspace(-0.8,1.1,50);
[P,~] = polyfit(clr,cdr,2);
cd_eval_root = polyval(P,cL_eval_root);


% NO VARIATION IN cd(y) (0012 the whole way)
for i = 1:length(alpha)
    % Finding cd at tip and root based on alpha 
    cdt = spline(cL_eval_tip,cd_eval_tip,(alpha(i)-2)/(2*pi));
    cdr = spline(cL_eval_tip,cd_eval_tip,alpha(i)/(2*pi));
    Cd0(i) = profile_drag_coeff(ct,cr,cdt,cdr,b,S,n);
    [~,c_L(i),c_Di(i)] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,(alpha(i) + 2)*pi/180,alpha(i)*pi/180,N);    
end

figure()
plot(alpha,Cd0,"Linewidth",1)
grid minor
xlabel("\alpha (deg)")
ylabel("C_{D,0}")
title("Profile Drag vs \alpha for whole aircraft")

% NO VARIATION IN cd(y) (0012 the whole way)
for i = 1:length(alpha)
    cdt = spline(cL_eval_tip,cd_eval_tip,(alpha(i)-2)/(2*pi));
    cdr = spline(cL_eval_root,cd_eval_root,alpha(i)/(2*pi));
    Cd0(i) = profile_drag_coeff(ct,cr,cdt,cdr,b,S,n);
    [~,c_L(i),c_Di(i)] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,(alpha(i) + 2)*pi/180,alpha(i)*pi/180,N);    
end

figure()
plot(alpha,Cd0,"Linewidth",1)
grid minor
xlabel("\alpha (deg)")
ylabel("C_{D,0}")
title("Profile Drag vs \alpha for whole aircraft")


%% CALCULATING TOTAL DRAG

% Finding drag
CD = Cd0 + c_Di;

% Plotting
figure()
hold on
plot(alpha,CD,"Linewidth",1)
plot(alpha,Cd0,"Linewidth",1)
plot(alpha,c_Di,"Linewidth",1)
grid minor
xlabel("\alpha (deg)")
ylabel("Coefficient of Drag (C_{D,0},C_{D,i},C_D)")
title("Total Drag vs \alpha")
legend("Total Drag","Profile Drag","Induced Drag","location","northwest")

close all
%% CALCULATING THRUST

% Finding density
h = 10000; %ft
h = h/3.281; %m
[~,~,~,rho] = atmosisa(h);
rho = rho/515.4*32.17; %slug/ft^3
W = 26000; %lbs

% Finding Velocity
V = sqrt(2*W./(c_L*rho*S));
V = V(V==real(V));

% Finding Drag
CD = CD(33:end);
T = CD.*(0.5*rho*V.^2*S);

% Calculating minimum thrust and associated V
[Tmin, Imin] = min(T)
Vmin = V(Imin)/1.688

%plotting results
figure()
hold on
plot(V/1.688,T,"Linewidth",1)
grid minor
xlabel("Airspeed (knots)")
ylabel("Thrust (lbf)")
xline(143.4)
legend("Thrust required","Max Speed","Location","Northwest")

figure()
hold on
plot(V/1.688,T,"Linewidth",1)
grid minor
xlabel("Airspeed (knots)")
ylabel("Thrust (lbf)")
xline(143.4)
xlim([0,145])
legend("Thrust required","Max Speed","Location","Northwest")
title("Thrust Required vs. Airspeed")
