
clear all
close all
clc
P = 1;
L = 1;
E = 1;
I = 1;

T = [6 -3*L; -3*L 5*L^2]

lambda = P*L^2/(60*E*I)
K_1 = [6 -3*L; -3*L 5*L^2];
K_G = [36 -3*L; -3*L 4*L^2];

[eigenVectors, eigenValues] = eig(K_1,K_G)
egenverdier = max(eigenValues)
Minste_egenverdi = min(egenverdier)
P_kr = Minste_egenverdi/ pi^2
P_euler = P_kr * lambda ;



