clear all
close all
clc

EI = 2.2*10^7;
L = 15;
m = 50;

x = linspace(0,L,1000);
modeShape = sin(pi*x/L);
modeShape_dd = -(pi/L)^2*sin(pi*x/L);

m_star = trapz(x, m*modeShape.^2); %integrerer 
k_star = trapz(x, EI*modeShape_dd.^2);

mu = 0.05;
omega_n = sqrt(k_star/m_star);

%Obtain the properties for the TMD

omega_n_TMD = omega_n * (1/(1+mu));
ksi_TMD = sqrt(3*mu/(8*(1+mu)));

m_TMD = m_star*mu;
k_TMD = m_TMD*omega_n_TMD^2;
c_TMD = ksi_TMD*2*m_TMD*omega_n_TMD;



