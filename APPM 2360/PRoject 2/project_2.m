% Author(s): Noah Vavoso, Andrew Patella, Karina Li
% Assignment title: APPM 2360 Project 2
% Creation date: 10/26/2023
%Purpose: Analyze probability and susceptibility to disease using Markov
%chains.

clear, clc, close all

%% Task A
Ps = 0.7;
Pe = 0.4;
Pi = 1;
Pr = 0.8;

%         S     E           I   R 
P_SEIR = [Ps,   Pe,         0,  1-Pr;
          1-Ps, 0,          0,  0;
          0,    0.5*(1-Pe), 0,  0;
          0,    0.5*(1-Pe), Pi, Pr];

% A3
S = [0; 1; 0; 0];
A3 = Markov(P_SEIR,S,1);

% A4
S = [1; 0; 0; 0];
A4 = Markov(P_SEIR,S,5);

%% Task B

S1 = [1; 0; 0; 0];
interval = 1:31;
PlotPls(P_SEIR,S1,interval)

% B2

S2=[.15;.85;0;0];
PlotPls(P_SEIR,S2,interval)

[eigenVectors,eigenValues] = eig(P_SEIR);





augmented_matrix = [eigenVectors S2];
sol_matrix = rref(augmented_matrix);
C = sol_matrix(:,5);

x_inf = C(3)*eigenVectors(:,3);

% B5
Pim = 0.5;
%.          S.    E.          I.     R.    Im
P_SEIRIm = [Ps,   Pe,         0,     1-Pr, 0;
            1-Ps, 0,          0,     0,    0;
            0,    0.5*(1-Pe), 0,     0,    0;
            0,    0.5*(1-Pe), 1-Pim, Pr,   0;
            0,    0,          Pim,   0,    1];

S3 = [1; 0; 0; 0; 0];
interval2 = 1:250;

PlotPls(P_SEIRIm,S3,interval2)


%% Task C
Pv = 0.25;
%.           S.    E.         I.     R.    Im.   V
P_SEIRVIm = [Ps,  Pe,         0,     1-Pr, 0,    0;
            1-Ps, 0,          0,     0,    0,    0;
            0,    0.5*(1-Pe), 0,     0,    0,    0;
            0,    0.5*(1-Pe), Pi,    Pr,   0,    0;
            0,    0,          0,     0,    1,    1-Pv;
            0,    0,          0,     0,    0,    Pv];
S4 = [.33; 0; 0; 0; 0; 0.67];

[c_eigVec,c_eigVal] = eig(P_SEIRVIm);


PlotPls(P_SEIRVIm,S4,interval)
    
  