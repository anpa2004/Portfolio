function PlotPls(P,S,interval)
    figure()
    lengthS = length(S);
    hold on;
    for i = interval
        MarkovStep = Markov(P,S,i)/norm(Markov(P,S,i),1);
        
        for j = 1:lengthS
            MarkovMatrix(j,i) = MarkovStep(j,1);
    
        end
    end
    
    for i = 1:lengthS
        plot(interval,MarkovMatrix(i,:),'LineWidth',.95)
    
    end
    grid on;
    legend('Suceptible','Exposed','Infected','Recovered','Immune','Vaccinated');
    xlabel('Time (days)');
    ylabel('Probability');
    title('Probability of Disease state');