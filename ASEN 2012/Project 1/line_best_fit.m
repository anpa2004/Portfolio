close all; clear; clc;

calibration_data = readmatrix('Static Test Stand Calibration Case 2.xlsx');

weight = calibration_data(:,1);
f_0_offset = calibration_data(:,2);
f_1_offset = calibration_data(:,3);
f_0_mv = calibration_data(:,4);
f_1_mv = calibration_data(:,5);


ftot_mv = f_0_mv+f_1_mv;

% f_long = cat(1,f_0_mv,f_1_mv);
% weight_long = cat(1,weight,weight);

f_adjusted_0 = weight.*(f_0_mv./(f_0_mv+f_1_mv));
f_adjusted_1 = weight.*(f_1_mv./(f_0_mv+f_1_mv));
f_adjusted_tot = f_adjusted_0+f_adjusted_1;

[p0,b0] = polyfit(f_adjusted_0,f_0_mv,1);
[p1,b1] = polyfit(f_adjusted_1,f_1_mv,1);
[ptot,btot] = polyfit(f_adjusted_tot,ftot_mv,1);

x = 0:.1:50;

[y0,delta0] = polyval(p0,x,b0);
[y1,delta1] = polyval(p0,x,b0);
[ytot,deltatot] = polyval(ptot,x,btot);

line0Fun = @(F) p0(1)*F + p0(2);
line1Fun = @(F) p1(1)*F + p1(2);

line0 = line0Fun(x);
line1 = line1Fun(x);



% figure(1);
% hold on;
% scatter(f_adjusted_1, f_0_mv);
% scatter(f_adjusted_0, f_1_mv);
% %plot(x,line0,x,line1,'blue');
% %plot(x,line0+2*delta0);
% plot(x,y0);
% plot(x,y1);
% plot(x,y0+2*delta0);
% plot(x,y0-2*delta1);
% legend('Data from Channel 0','Data from Channel 1','Best Fit 0','Best fit 1','','Location','Southeast');
% 
% grid on;
% xlabel('Adjusted Force (lbf)');
% ylabel('Measured Load (mV)');
% hold off;

figure(2);
hold on;
scatter(f_adjusted_tot,ftot_mv);
plot(x,ytot,'color','blue')
plot(x,ytot+2*deltatot,'color','#D95319');
plot(x,ytot-2*deltatot,'color',"#D95319");
grid on;
xlabel('Adjusted Weight (lbf)');
ylabel('Measured Load (mV)');
title('Adjusted mV Output for Applied Load');
legend('Calibration data','Best Fit line','F best \pm 2\sigma','location','Northwest')

values = cat(1,p0,p1);
writematrix(values,'linear_fit_parameters.csv');

% figure(2);
% hold on;
% scatter(weight_long,f_long);
% grid on;

% for i =1:length(weight)
%     fprintf('Weight = %f \n',weight(i));
%     fprintf('Sum of adj = %f \n',f_adjusted_1(i)+f_adjusted_0(i));
% end

