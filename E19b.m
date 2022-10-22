clear all
close all
clc

K = [4.92 -4.8 0 0 0; -4.8 32 8 0 0; 0 8 48 -4.8 8; 0 0 -4.8 1.92 -4.8; 0 0 8 -4.8 16]
M = [21.5714 -3.9286 0 0 0; -3.9286 7.1429 -2.6786 0 0; 
    0 -2.6786 10.7143 2.3214 -2.6786; 0 0 2.3214 5.5714 -3.9286; 0 0 -2.6786 -3.9286 3.5714]

[eigenVectors, eigenValues] = eig(K,M)
omega_n = sqrt(diag(eigenValues))

numDofs = size(K);
temporary_mass_matrix = transpose(eigenVectors)*M*eigenVectors;
massNormalized_eigenVectors = zeros(numDofs);
for i = 1:numDofs(1)
    massNormalized_eigenVectors(:,i) = eigenVectors(:,i)/temporary_mass_matrix(i,i);
end

%calculate the displacement normalized eigenvectors
dispNormalized_eigVecs = zeros(numDofs);
for i = 1:numDofs(1)
    dispNormalized_eigVecs(:,i) = eigenVectors(:,i)/eigenVectors(i,i);
end
Phi = eigenVectors

K_stjerne = transpose(Phi)*K*Phi
M_stjerne = transpose(Phi)*M*Phi

xi_1 = 0.03;
xi_3 = 0.05;

omega_1 = omega_n(1);
omega_3 = omega_n(3);

Dampingratiovector = [xi_1; xi_3];
Naturalfrequency = 0.5* [1/omega_1 omega_1; 1/omega_3 omega_3];
koeffisienter = inv(Naturalfrequency) * Dampingratiovector
C_stjerne = koeffisienter(1)*M_stjerne + koeffisienter(2)*K_stjerne


P_0 = [1; 1; 1; 1; 1];
P_0_stjerne = transpose(Phi)*P_0


u_0 = [1; 0; 0; 0; 0];
u_0_dot = [0; 0; 0; 0; 0];

y_0 = inv(Phi)*u_0;
y_0_dot = inv(Phi)*u_0_dot;

numDofs = size(K_stjerne);

%Obtain the damping ratio in each mode from the Rayleigh damping:

xis = zeros(1,numDofs(1));
for i = 1:numDofs(1)
    xis(i) = C_stjerne(i,i) / (2*M_stjerne(i,i)* omega_n(i));
end

%Calculate the contribution from all modes:
% Both h and p solution
omega_0 = 0.32 

timelength = 10000;
t = linspace(0,2*pi/omega_0*35, timelength);
%create a matrix of correct shape to store our modal solutions:

y = zeros([numDofs(1), timelength]);
for i=1:numDofs
    beta =  omega_0 / omega_n(i);
    D = 1/sqrt((1-beta^2)^2 + (2*xis(i)*beta)^2);
    theta = atan2(2*xis(i)*beta, 1-beta^2);
    y_p = P_0_stjerne(i)/K_stjerne(i,i)*D*cos(omega_0*t - theta);
    omega_D = omega_n(i)*sqrt(1-xis(i)^2);
    C_1 = y_0(i) - P_0_stjerne(i)/K_stjerne(i,i)*D*cos(-theta);
    C_2 = P_0_stjerne(i)/K_stjerne(i,i)*D*sin(-theta) + y_0_dot(i) + C_1*xis(i)*omega_n(i);
    y_h = exp(-xis(i)*omega_n(i)*t) .* (C_1*cos(omega_D*t) + C_2*sin(omega_D*t));
    
    y_mode = y_h + y_p;
    y(i,:) = y_mode;
end

r = Phi*y;
plot(t,r(1,:));

















