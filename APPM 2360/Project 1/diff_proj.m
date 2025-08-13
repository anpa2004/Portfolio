close all; clear; clc;

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