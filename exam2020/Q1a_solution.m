close all
clear all
clc;
%% Parameters
EI = 100;       %Bending Stifness
L = 2;          %Length
m_c =  2;       %Column Mass
m_s = 10;       %Floor Mass

%% Stiffness and Mass Matrix
% K-matrisen er funnet ved bruk av formell metode, K = transpose(a)*k*a
a_AB = [1 0 0 0; 
        0 0 0 0; 
        0 0 0 0; 
        0 0 0 0];
a_BC = [0 1 0 0;
        0 0 0 0; 
        1 0 0 0; 
        0 0 0 0];
a_CD = [0 0 1 0; 
        0 0 0 0; 
        0 1 0 0; 
        0 0 0 0];

% Geometri
L_AB = L; 
L_BC = L;
L_CD = L;
EI_AB = 2*EI; % to søyler
EI_BC = 2*EI;
EI_CD = 2*EI;

k_AB = EI_AB/L_AB^3 * [12 -6*L_AB -12 -6*L_AB; -6*L_AB 4*L_AB^2 6*L_AB 2*L_AB^2; -12 6*L_AB 12 6*L_AB;-6*L_AB 2*L_AB^2 6*L_AB 4*L_AB^2];
k_BC = EI_BC/L_BC^3 * [12 -6*L_BC -12 -6*L_BC; -6*L_BC 4*L_BC^2 6*L_BC 2*L_BC^2; -12 6*L_BC 12 6*L_BC;-6*L_BC 2*L_BC^2 6*L_BC 4*L_BC^2];
k_CD = EI_CD/L_CD^3 * [12 -6*L_CD -12 -6*L_CD; -6*L_CD 4*L_CD^2 6*L_CD 2*L_CD^2; -12 6*L_CD 12 6*L_CD;-6*L_CD 2*L_CD^2 6*L_CD 4*L_CD^2];

K = (transpose(a_AB)*k_AB*a_AB + transpose(a_BC)*k_BC*a_BC + transpose(a_CD)*k_CD*a_CD) ;
K(:,4) = [];
K(4,:) = []

% Mass matrix:
M1 = m_c*L/420*[156 0 0; 0 0 0; 0 0 0];
M2 = m_c*L/420*[156 54 0; 54 156 0; 0 0 0];
M3 = m_c*L/420*[0 0 0; 0 156 54; 0 54 156];


%M2 = m_c*L/420*[156
%M3 = m_c*L/420*[156

%Mass matrix for the floor by unit acceleration:
m_floor = m_s*L/420*[420 0 0; 0 420 0; 0 0 420]; % 

M = 2*M1 + 2*M2 + 2*M3 + m_floor

