%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrix_project.m is designed to take in a shape of a tmv graph and
% condition it to reasonable values. Then it calculates APD vales and
% outputs the conditioned data and APD values to be used in later
% functions. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc;

%Reading in data
data = readmatrix("APD_graph.xlsx");

time = data(:,1);
v = data(:,2);

%Scaling and conditioning to match image (Approximate)
time = time*290/(max(time));
v = v*75/abs(max(v));
v = v-85;

%% Calculating APD values
first = 4;
last = 9;

% APD^[] values that you wish to measure
values = first:1:last;

%The change in voltage across a heartbeat
v_diff = abs(max(v)-min(v));

%Splitting data into rise and fall for ease
[peakV,peaki] = max(v);
riseV = v(1:peaki);
riseT = time(1:peaki);

fallV = v(peaki+9:end);
fallT = time(peaki+9:end);

t_stim = time(1); 

% Calculating different APD values
for i = 1:length(values)
    v_apd{values(i)} = max(v)-values(i)/10*v_diff;

    t1{values(i)} = interp1(riseV,riseT,v_apd{values(i)});
    t2{values(i)} = interp1(fallV,fallT,v_apd{values(i)});

    APDv{values(i)} = t2{values(i)}-t_stim;
end


for i = first:last
    APD(i-(first-1)) = APDv{i};
end

csvwrite("tmv.csv",v)
csvwrite("time.csv",time)
csvwrite("APD_ss",APD);

% Plots
figure()
hold on;
grid minor;
plot(time,v,"Linewidth",2);
xlabel("Time (ms)");
ylabel("Transmembrane Voltage (mV)");

for i = first:last
    yline(v_apd{i});
    
end

xline(t1{i},"Color",'r');
text(45,v_apd{4}(1,1)+2.5,"APD^{40}");
text(55,v_apd{5}(1,1)+2.5,"APD^{50}");
text(65,v_apd{6}(1,1)+2.5,"APD^{60}");
text(75,v_apd{7}(1,1)+2.5,"APD^{70}");
text(85,v_apd{8}(1,1)+2.5,"APD^{80}");
text(95,v_apd{9}(1,1)+2.5,"APD^{90}");
title("Transmembrane Voltage for Heart (Model Data)")
legend('Transmembrane Voltage Data','APD^{[ ]} values','','','','','','Stimulus')


figure()
plot(10*values,APD,'linewidth',1);
grid minor;
xlabel('APD^{[]}');
ylabel('APD [ms]')
title("APD values for Different Percentages of Restabilization")


