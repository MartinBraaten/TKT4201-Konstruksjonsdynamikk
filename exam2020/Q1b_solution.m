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

% Matlab makes the eigenvectors mass normalized
[eigenVectors, eigenValues] = eig(K,M)
omega_n = sqrt(diag(eigenValues)) %Natural frequencies

mode_shapes = eigenVectors
Phi = eigenVectors;
Phi(:,1) = Phi(:,1)/Phi(1,1);
Phi(:,2) = Phi(:,2)/Phi(1,2);
Phi(:,3) = Phi(:,3)/Phi(1,3)
% Matlab makes the eigenvectors mass normalized
% To check if this is true, tranpose(eigenVectors)*M*eigenVectors == identity matrix
eigenVectors'*M*eigenVectors





