clear; close all; clc;

%% Section C
%Variable values (initial and final values)
t_0 = 0;
t_f = 24;
h = 0.1;
m_0 = 75;
eT_0=25;
ek=0.25;

%Declaration of t array for function
et=zeros((t_f-t_0)/h,1);   
etime = t_0:h:t_f;

%Function handle for the differential equation
ef = @(et,eT).25*(m_0-eT);

%Using the Runga Kutta function
[eint_t, eint_T] = rk4(t_0,t_f,240,eT_0,ef);

%Calculated result
T_of_t = m_0+(eT_0-m_0)*exp(-ek*etime);

%Plotting results
figure (1);
subplot(1,2,1);
hold on;
plot(eint_t,eint_T,'color','r','LineStyle','--');
plot(etime,T_of_t,'color','b');
legend('Runga Kutta Aproximation','Exact Solution');
grid on;
xlabel('Time');
ylabel('Temperature (Ferenheit)');
title('Temperature change from outside temperature (M_0)')
hold off;

subplot(1,2,2);
hold on;
plot(etime,abs(T_of_t - eint_T));
title('Error of Runga Kutta Aproximation');
xlabel('Time');
ylabel('Absolute Error (^oF)');
grid on;

%% Section D

% plot D.1================================================================
dx=1:.1:24;
dy=(28/3)*(atan(sinh(0.75*dx-7.5)))+79.65044;

figure(2)
subplot(1,2,1)
plot(dx,dy,'-b')
grid on
title('Building Temperature Due to People, Lights, and Machinery')
xlabel('Time (hrs)')
ylabel('Temperature (^oF)')

% plot D.2
dy2=7*(sech(0.75*dx-7.5));

subplot(1,2,2)
plot(dx,dy2,'-b')
grid on
title('Rate of Change of Temperature Due to People, Lights, and Machinery')
xlabel('Time (hrs)')
ylabel('Rate of Change of Temperature (^oF/hrs)')

% D.3
[dm,dt]=max(dy);
dmin=(dt/10-floor(dt/10))*60;
dtime=floor(dt/10);

fprintf('The maximum temperature in the building is %3.2f F and occurs at %g hours %g minutes.\n',dm,dtime,dmin);

%========================================================
%% Section F

%Variable values (initial and final values)--------------------------------
t_0 = 0;
t_f = 24;
h = 0.1;
T_0=75;
K = 0.5;
t_f2 = 72.0;
time = t_0:h:t_f2;

%Declaration of t array for function---------------------------------------
t=zeros((t_f-t_0)/h,1);   

%Function handle for the differential equation-----------------------------
f = @(t,T)7*sech(.75*(t-10))+2*(77-T);
g = @(t,T)0.25*(85-10*cos(pi*(t-5)/12)-T);
j = @(t,T)0.25*(85-10*cos(pi*(t-5)/12)-T)+K*(77-T);
l = @(t,T)0.25*(85-10*cos(pi*(t-5)/12)-T)+K*(77-T)+7*sech(.75*(t-10))+2*(77-T);
M = @(t,T)85-10*cos(pi*(t-5)/12);

%Using the Runga Kutta function--------------------------------------------
[int_t, int_T] = rk4(t_0,t_f,240,T_0,f);
[int_t1, int_T1] = rk4(t_0,t_f,240,T_0,g);
[int_t2, int_T2] = rk4(t_0,t_f,240,T_0,j);
[int_t3, int_T3] = rk4(t_0,t_f2,240,T_0,l);

%plotting all the values---------------------------------------------------
figure(3)
subplot(2,2,1);
plot(int_t,int_T);
title('Aproximation of temperature (a)');
xlabel('Time (hours)');
ylabel('Temperature (^oF)');
yline(81,'r--');
legend('Building Temperature','Danger Temperature');
grid on;

subplot(2,2,2);
plot(int_t1,int_T1);
title('Aproximation of temperature (b)');
xlabel('Time (hours)');
ylabel('Temperature (^oF)');
yline(81,'r--');
legend('Building Temperature','Danger Temperature');
grid on;

subplot(2,2,3);
plot(int_t2,int_T2);
title('Aproximation of temperature (c)');
xlabel('Time (hours)');
ylabel('Temperature (^oF)');
yline(81,'r--');
legend('Building Temperature','Danger Temperature');
grid on;

subplot(2,2,4);
hold on;
plot(int_t3,int_T3,time,M(time));
title('Aproximation of temperature (d)');
xlabel('Time (hours)');
ylabel('Temperature (^oF)');
yline(81,'r--');
legend('Building Temperature','Outside Temperature','Danger Temperature');
grid on;

%Calculating max values----------------------------------------------------
[max_T,max_t] = max(int_T);
[max_T1,max_t1] = max(int_T1);
[max_T2,max_t2] = max(int_T2);
[max_T3,max_t3] = max(int_T3);

%reducing to hour/mins-----------------------------------------------------
actual1 = (max_t/10);
rem = actual1 - floor(actual1);
actual1m = rem*60;

actual2 = max_t1/10;
rem = actual2 - floor(actual2);
actual2m = rem*60;

actual3 = max_t2/10;
rem = actual3 - floor(actual3);
actual3m = rem*60;

actual4 = max_t3/10;
rem = actual4 - floor(actual4);
actual4m = rem*60;

%Displaying Max values-----------------------------------------------------
fprintf('Maximum Value 1: %2.2f F at %g hours %g minutes\n',max(int_T),floor(actual1),actual1m);
fprintf('Maximum Value 2: %2.2f F at %g hours %g minutes \n',max(int_T1),floor(actual2),actual2m);
fprintf('Maximum Value 3: %2.2f F at %g hours %g minutes\n',max(int_T2),floor(actual3),actual3m);
fprintf('Maximum Value 4: %2.2f F at %g hours %g minutes\n',max(int_T3),floor(actual4),actual4m);




