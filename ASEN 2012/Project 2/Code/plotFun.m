function plotFun(state,parameter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function takes in cells of time and state from ode45 and plots all
% of the struct rows on similar graphs

% No return is necessary, plots will just be created


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot of trajectory 
    figure()
    hold on;
    
    %For loop to plot all of the 8 different lines
    for i = 1:8
        plot(state{i}(:,1),state{i}(:,3),'Linewidth',1);
    end

    grid on;
    xlabel ('Distance (m) ');
    ylabel ( 'Distance (m)');


    text = parameter;

    title(sprintf('Trajectory based on changing %s',text));
    legend(sprintf('%s 1',text),sprintf('%s 2',text), ...
        sprintf('%s 3',text),sprintf('%s 4',text), ...
        sprintf('%s 5',text),sprintf('%s 6',text), ...
        sprintf('%s 7',text),sprintf('%s 8',text));

    hold off;

end