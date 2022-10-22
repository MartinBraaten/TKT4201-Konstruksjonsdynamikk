%Script to plot the response of the two physical degrees of freedom
clear all
close all
clc

dt = 0.01; % time decrement
t = 0:dt:100; % time interval
t = [0:0.01:100];

%f = 0.8337;  % low
% f = 1.562;  % between
 f = 2.681;  % high
% f = 0.9295298;  % first resonance
% f = 1.90384326;  % second resonance

w_load = f*2*pi;   % brukt i partikulærløsning
P0 = 0.3236*w_load^2; 
P = [-0.0423*P0, 0.0448*P0];  % P* = phi*P0
w = [5.8404, 11.9622];    % brukt i komplementærløsning
Beta = [w_load/5.8404, w_load/11.9622];
K = [34.1105, 143.0932]; % w^2
xi = [0.002671, 0.001149]; % Funnet fra decay_hor_deltaM og Decay_rot

D1 = P(1)/K(1) *1/((1-Beta(1)^2)^2 + (2*xi(1)*Beta(1))^2); % D er nesten det samme som G
D2 = P(2)/K(2) *1/((1-Beta(2)^2)^2 + (2*xi(2)*Beta(2))^2);

A1 = D1 *2*xi(1)*Beta(1); % A og B er konstanter i komplementær løsning
A2 = D2 *2*xi(2)*Beta(2);
B1 = D1 *1/sqrt(1-xi(1)^2) *(2*xi(1)^2*Beta(1) - Beta(1)*(1-Beta(1)^2));
B2 = D2 *1/sqrt(1-xi(2)^2) *(2*xi(2)^2*Beta(2) - Beta(2)*(1-Beta(2)^2));

y_p1 = D1* ((1-Beta(1)^2) *sin(w_load*t) - 2*xi(1)*Beta(1) *cos(w_load*t)); % Partikulær løsning
y_p2 = D2* ((1-Beta(2)^2) *sin(w_load*t) - 2*xi(2)*Beta(2) *cos(w_load*t));

y_c1 = exp(-xi(1)*w(1)*t) .* (A1 *cos(w(1)*sqrt(1-xi(1)^2)*t)) + (B1 *sin(w(1)*sqrt(1-xi(1)^2)*t)); % Komplementær
y_c2 = exp(-xi(2)*w(2)*t) .* (A2 *cos(w(2)*sqrt(1-xi(2)^2)*t)) + (B2 *sin(w(2)*sqrt(1-xi(2)^2)*t));

y1 = y_p1 + y_c1;
y2 = y_p2 + y_c2;

u1 = -0.0401*y1 - 0.0067*y2;    % Phi*y
u2 = -0.0014*y1 + 0.0343*y2;


% plot
figure();
subplot(2,1,1);
plot(t, u1*1000);
title('transversal response of midpoint');
xlabel('seconds');
ylabel('mm');
subplot(2,1,2);
plot(t, u2);
title('rotational response')
xlabel('seconds');
ylabel('radians');
sgtitle('Load frequency: 2.681 Hz');