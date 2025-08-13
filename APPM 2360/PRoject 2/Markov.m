function [probability] = Markov(p,s,t)
    probability = (p^t)*s;
end