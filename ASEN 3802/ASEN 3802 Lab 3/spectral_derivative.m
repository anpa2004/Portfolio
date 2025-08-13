function ddx = spectral_derivative(x,y,n)
% This function takes in a data set y=y(x) and estimates the nth
% derivative of the dataset using spectral methods
    L = 2*x(end);
    n = length(x);
    kappa = (2*pi/L)*[-n/2:n/2-1];
    kappa = fftshift(kappa);
    yhat = fft(y);
    ddxhat = yhat.*(1i*kappa).^n;
    ddx = ifft(ddxhat);