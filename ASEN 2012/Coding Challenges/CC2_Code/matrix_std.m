function [wMean,wStd] = matrix_std(A)
sigma = std(A);

xbar = mean(A);

W = 1./sigma.^2;

wMean = sum(W.*xbar)./sum(W);

wStd = 1./sqrt(sum(W));

end