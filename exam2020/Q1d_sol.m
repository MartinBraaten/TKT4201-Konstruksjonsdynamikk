close all
clear all
clc;
%% Parameters
EI = 100;       %Bending Stifness
L = 2;          %Length
m_c =  2;       %Column Mass
m_s = 10;       %Floor Mass

%% Stiffness and Mass Matrix
K = [800  -400     0
    -400   800  -400
       0  -400   450];
M = [25    3    0
      3   25    3
      0    3   25];

%%
% Natural Frequency and shape modes
w_n =[0.2816
      0.8187
      1.2252];
modes=[   0.0322    0.0712   -0.0525
          0.0570    0.0180    0.0733
          0.0685   -0.0667   -0.0499];
      
%% Forced Vibration Response
p_0=100;        %Force Amplitude

w = linspace(0,3, 10000)
t = linspace(0,10, 10000)
P = p_0 * sin(w*t)*[1, 0, 0] % Harmonic load

numDofs = size(K);
K_star = modes'*K*modes
M_star = modes'*M*modes
P_star = modes'*P
for i = 1:numDofs:
    beta = w/w_n(i)
    theta = atan2(1/(1-beta^2))
    D = 1/sqrt((1-beta^2)^2)







