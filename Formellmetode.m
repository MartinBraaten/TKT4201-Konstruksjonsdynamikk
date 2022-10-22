clear all
close all
clc

syms P L EI

a_AB = [0 0 0 0; 0 0 0 0; 1 0 0 0; 0 1 0 0]
a_BC = [0 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 1 0]
a_CD = [-1 0 0 0; 0 0 1 0; 0 0 0 0; -3/(2*L) 0 -0.5 0]

%Geometri
L_AB = L; 
L_BC = L;
L_CD = L;
EI_AB = 2*EI; 
EI_BC = EI;
EI_CD = EI;
N_AB = 0 ; 
N_BC = -P ; 
N_CD = 0 ;
Factor_K_1 = EI/L^3 ;
Factor_K_G = -P/(30*L);

% Er noen av disse = 0?
k_G_AB = 0 ;%N_AB/(30*L_AB) * [36 -3*L_AB -36 -3*L_AB; -3*L_AB 4*L_AB^2 3*L_AB -L_AB^2; -36 3*L_AB 36 3*L_AB; -3*L_AB -L_AB^2 3*L_AB 4*L_AB^2];
k_G_BC = N_BC/(30*L_BC) * [36 -3*L_BC -36 -3*L_BC; -3*L_BC 4*L_BC^2 3*L_BC -L_BC^2; -36 3*L_BC 36 3*L_BC; -3*L_BC -L_BC^2 3*L_BC 4*L_BC^2];
k_G_CD = 0 ;%N_CD/(30*L_CD) * [36 -3*L_CD -36 -3*L_CD; -3*L_CD 4*L_CD^2 3*L_CD -L_CD^2; -36 3*L_CD 36 3*L_CD; -3*L_CD -L_CD^2 3*L_CD 4*L_CD^2];

% 3Dof bjelkeelement
% k_AB = 3*EI_AB/L_AB^3 * [1 -L_AB -1 0; -L_AB L_AB^2 L_AB 0; -1 L_AB  0;0 0 0 0];
% k_BC = 3*EI_BC/L_BC^3 * [1 -L_BC -1 0; -L_BC L_BC^2 L_BC 0; - L_BC 1 0;0 0 0 0];
k_CD = 3*EI_CD/L_CD^3 * [1 -L_CD -1 0; -L_CD L_CD^2 L_CD 0; -1 L_CD 1 0;0 0 0 0]

% 4Dof bjelkeelement
k_AB = EI_AB/L_AB^3 * [12 -6*L_AB -12 -6*L_AB; -6*L_AB 4*L_AB^2 6*L_AB 2*L_AB^2; -12 6*L_AB 12 6*L_AB;-6*L_AB 2*L_AB^2 6*L_AB 4*L_AB^2];
k_BC = EI_BC/L_BC^3 * [12 -6*L_BC -12 -6*L_BC; -6*L_BC 4*L_BC^2 6*L_BC 2*L_BC^2; -12 6*L_BC 12 6*L_BC;-6*L_BC 2*L_BC^2 6*L_BC 4*L_BC^2];
% k_CD = EI_CD/L_CD^3 * [12 -6*L_CD -12 -6*L_CD; -6*L_CD 4*L_CD^2 6*L_CD 2*L_CD^2; -12 6*L_CD 12 6*L_CD;-6*L_CD 2*L_CD^2 6*L_CD 4*L_CD^2];

K_1 = transpose(a_AB)*k_AB*a_AB + transpose(a_BC)*k_BC*a_BC + transpose(a_CD)*k_CD*a_CD;
K_G = transpose(a_AB)*k_G_AB*a_AB + transpose(a_BC)*k_G_BC*a_BC + transpose(a_CD)*k_G_CD*a_CD;
K_2 = K_1 + K_G ;

K_1_cor = K_1/Factor_K_1
K_G_cor = K_G/Factor_K_G;

%jevnt fordelt belastning
p_AB = 0 ;
p_BC = 0;
p_CD = 0 ;
%Kraft i midtspenn
P_AB = 0;
P_BC = 0;
P_CD = -P;
%R_0_faktor = P*L/16

% 3Dof bjelkeelement
% S_0_AB = p_AB*L_AB/8*[-5; L_AB; -3; 0] + P_AB/16*[-11; 3*L_AB; -5; 0];
% S_0_BC = p_BC*L_BC/8*[-5; L_BC; -3; 0] + P_BC/16*[-11; 3*L_BC; -5; 0];
S_0_CD = p_CD*L_CD/8*[-5; L_CD; -3; 0] + P_CD/16*[-11; 3*L_CD; -5; 0]
% 4Dof bjelkeelement
S_0_AB = p_AB*L_AB/12*[-6; L_AB; -6; -L_AB] + P_AB/8*[-4; L_AB; -4; -L_AB];
S_0_BC = p_BC*L_BC/12*[-6; L_BC; -6; -L_BC] + P_BC/8*[-4; L_BC; -4; -L_BC];
% S_0_CD = p_CD*L_CD/12*[-6; L_CD; -6; -L_CD] + P_CD/8*[-4; L_CD; -4; -L_CD];

R_0_AB = transpose(a_AB)*S_0_AB;
R_0_BC = transpose(a_BC)*S_0_BC;
R_0_CD = transpose(a_CD)*S_0_CD

R_0 = (R_0_AB + R_0_BC + R_0_CD) %/R_0_faktor
R_K = [-P; 0; 0; 0];
R = R_K - R_0

S_faktor = P/44 ;
r = P*L^2/(528*EI)*[-59*L; 60; -6; 0];

S_AB = ((k_AB*a_AB*r) + (S_0_AB))/S_faktor
S_BC = ((k_BC*a_BC*r) + (S_0_BC))/S_faktor
S_CD = ((k_CD*a_CD*r) + (S_0_CD))/S_faktor


% Knekningsmatrise
T_A = [1 0; 0 1; 0 1];
T_Sym = [0; 1; -1];
K_1_A = transpose(T_A)*K_1_cor*T_A;
K_G_A = transpose(T_A)*K_G_cor*T_A;
K_1_Sym = transpose(T_Sym)*K_1_cor*T_Sym;
K_G_Sym = transpose(T_Sym)*K_G_cor*T_Sym;

disp('husk å sette inn disse faktorene foran matrisen')
% Faktorisere K_1_A
f_1 = factor(abs(K_1_A(1)));
f_2 = factor(abs(K_1_A(2)/L));
f_3 = factor(abs(K_1_A(3)/L));
f_4 = factor(abs(K_1_A(4)/L^2));
i_1 = intersect(f_1,f_2);
i_2 = intersect(f_3,f_4);
i_3 = intersect(i_1,i_2)
if isempty(i_3) 
    i_3 = 1;
end

% Faktorisere K_1_G
r_1 = factor(abs(K_G_A(1)));
r_2 = factor(abs(K_G_A(2)/L));
r_3 = factor(abs(K_G_A(3)/L));
r_4 = factor(abs(K_G_A(4)/L^2));
j_1 = intersect(r_1,r_2);
j_2 = intersect(r_3,r_4);
j_3 = intersect(j_1,j_2)
if isempty(j_3) 
    j_3 = 1;
end

K_1_A_factorized = K_1_A / i_3
K_G_A_factorized = K_G_A / j_3

Lambda = -Factor_K_G*j_3 / (Factor_K_1*i_3)
KG = double(K_G_A_factorized.*[1 1/L; 1/L 1/L^2])
K1 = double(K_1_A_factorized.*[1 1/L; 1/L 1/L^2])

[eigenVectors, eigenValues] = eig(K1,KG)
egenverdier = max(eigenValues)
Minste_egenverdi = min(egenverdier)
P_kr = Minste_egenverdi/ pi^2
digits(5);
P_euler = vpa(P_kr / Lambda) 

if abs(eigenValues(1)) < abs(eigenValues(4))
    eigenVectors(1)/eigenVectors(2)
else 
    eigenVectors(4)/eigenVectors(3)
end


