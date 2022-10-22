clear all
close all
clc
syms P L EI p c s

%theta = deg2rad(3/4) ; %Fra node 1 til 2 fra global x-akse
c = 2/sqrt(5);
s = -1/sqrt(5);
T = [c^2 c*s -c^2 -c*s; c*s s^2 -c*s -s^2; -c^2 -c*s c^2 c*s; -c*s -s^2 c*s s^2];

% [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]
a_AB = [0 0 0 0; 
        0 0 0 0; 
        0 0 0 0; 
        0 0 0 0]
a_BC = [0 0 0 0;
        0 0 0 0; 
        0 0 0 0; 
        0 0 0 0]
a_CD = [0 0 0 0; 
        0 0 0 0; 
        0 0 0 0; 
        0 0 0 0]
a_DE = [0 0 0 0; 
        0 0 0 0; 
        0 0 0 0; 
        0 0 0 0]
%Geometri
L_AB = L; 
L_BC = L;
L_CD = L;
L_DE = L;
EI_AB = EI; 
EI_BC = EI;
EI_CD = EI;
EI_DE = EI;
%Normalkraft i søylen
N_AB = 0; 
N_BC = 0; 
N_CD = 0;
N_DE = 0;
%jevnt fordelt belastning
p_AB = 0;
p_BC = 0;
p_CD = 0;
p_DE = 0;
%Kraft i midtspenn
P_AB = 0;
P_BC = 0;
P_CD = 0;
P_DE = 0;
R_0_faktor = 1 %P*L/1
Factor_K_G = -P/(L);
Factor_K_1 = EI/(2*L^3);
R_K = [0; 0; 0; 0]; % [0; 0; 0; 0]
r = 1/(EI)*[0; 0; 0; 0] % [0; 0; 0; 0]


% 3Dof bjelkeelement
 k_AB = 3*EI_AB/L_AB^3 * [1 -L_AB -1 0; -L_AB L_AB^2 L_AB 0; -1 L_AB 1 0; 0 0 0 0];
 k_BC = 3*EI_BC/L_BC^3 * [1 -L_BC -1 0; -L_BC L_BC^2 L_BC 0; -1 L_BC 1 0; 0 0 0 0];
% k_CD = 3*EI_CD/L_CD^3 * [1 -L_CD -1 0; -L_CD L_CD^2 L_CD 0; -1 L_CD 1 0; 0 0 0 0];
 k_DE = 3*EI_DE/L_DE^3 * [1 -L_DE -1 0; -L_DE L_DE^2 L_DE 0; -1 L_DE 1 0; 0 0 0 0];

% 4Dof bjelkeelement
%k_AB = EI_AB/L_AB^3 * [12 -6*L_AB -12 -6*L_AB; -6*L_AB 4*L_AB^2 6*L_AB 2*L_AB^2; -12 6*L_AB 12 6*L_AB;-6*L_AB 2*L_AB^2 6*L_AB 4*L_AB^2];
%k_BC = EI_BC/L_BC^3 * [12 -6*L_BC -12 -6*L_BC; -6*L_BC 4*L_BC^2 6*L_BC 2*L_BC^2; -12 6*L_BC 12 6*L_BC;-6*L_BC 2*L_BC^2 6*L_BC 4*L_BC^2];
%k_CD = EI_CD/L_CD^3 * [12 -6*L_CD -12 -6*L_CD; -6*L_CD 4*L_CD^2 6*L_CD 2*L_CD^2; -12 6*L_CD 12 6*L_CD;-6*L_CD 2*L_CD^2 6*L_CD 4*L_CD^2];
%k_DE = EI_DE/L_DE^3 * [12 -6*L_DE -12 -6*L_DE; -6*L_DE 4*L_DE^2 6*L_DE 2*L_DE^2; -12 6*L_DE 12 6*L_DE;-6*L_DE 2*L_DE^2 6*L_DE 4*L_DE^2];

k_G_AB = N_AB/(30*L_AB) * [36 -3*L_AB -36 -3*L_AB; -3*L_AB 4*L_AB^2 3*L_AB -L_AB^2; -36 3*L_AB 36 3*L_AB; -3*L_AB -L_AB^2 3*L_AB 4*L_AB^2];
k_G_BC = N_BC/(30*L_BC) * [36 -3*L_BC -36 -3*L_BC; -3*L_BC 4*L_BC^2 3*L_BC -L_BC^2; -36 3*L_BC 36 3*L_BC; -3*L_BC -L_BC^2 3*L_BC 4*L_BC^2];
k_G_CD = N_CD/(30*L_CD) * [36 -3*L_CD -36 -3*L_CD; -3*L_CD 4*L_CD^2 3*L_CD -L_CD^2; -36 3*L_CD 36 3*L_CD; -3*L_CD -L_CD^2 3*L_CD 4*L_CD^2];
k_G_DE = N_DE/(30*L_DE) * [36 -3*L_DE -36 -3*L_DE; -3*L_DE 4*L_DE^2 3*L_DE -L_DE^2; -36 3*L_DE 36 3*L_DE; -3*L_DE -L_DE^2 3*L_DE 4*L_DE^2];

K_1 = transpose(a_AB)*k_AB*a_AB + transpose(a_BC)*k_BC*a_BC + transpose(a_CD)*k_CD*a_CD + transpose(a_DE)*k_DE*a_DE
K_G = transpose(a_AB)*k_G_AB*a_AB + transpose(a_BC)*k_G_BC*a_BC + transpose(a_CD)*k_G_CD*a_CD + transpose(a_DE)*k_G_DE*a_DE;
K_1_cor = K_1/Factor_K_1
K_G_cor = K_G/Factor_K_G


% 3Dof bjelkeelement
%S_0_AB = p_AB*L_AB/8*[-5; L_AB; -3; 0] + P_AB/16*[-11; 3*L_AB; -5; 0];
%S_0_BC = p_BC*L_BC/8*[-5; L_BC; -3; 0] + P_BC/16*[-11; 3*L_BC; -5; 0];
%S_0_CD = p_CD*L_CD/8*[-5; L_CD; -3; 0] + P_CD/16*[-11; 3*L_CD; -5; 0];
%S_0_DE = p_DE*L_DE/8*[-5; L_DE; -3; 0] + P_DE/16*[-11; 3*L_DE; -5; 0];

% 4Dof bjelkeelement
S_0_AB = p_AB*L_AB/12*[-6; L_AB; -6; -L_AB] + P_AB/8*[-4; L_AB; -4; -L_AB];
S_0_BC = p_BC*L_BC/12*[-6; L_BC; -6; -L_BC] + P_BC/8*[-4; L_BC; -4; -L_BC];
S_0_CD = p_CD*L_CD/12*[-6; L_CD; -6; -L_CD] + P_CD/8*[-4; L_CD; -4; -L_CD];
S_0_DE = p_DE*L_DE/12*[-6; L_DE; -6; -L_DE] + P_DE/8*[-4; L_DE; -4; -L_DE];

R_0_AB = transpose(a_AB)*S_0_AB
R_0_BC = transpose(a_BC)*S_0_BC
R_0_CD = transpose(a_CD)*S_0_CD
R_0_DE = transpose(a_DE)*S_0_DE;

R_0 = (R_0_AB + R_0_BC + R_0_CD + R_0_DE) %/R_0_faktor
R = R_K - R_0

S_faktor = 1 ;

S_AB = ((k_AB*a_AB*r) + (S_0_AB))/S_faktor
S_BC = ((k_BC*a_BC*r) + (S_0_BC))/S_faktor
S_CD = ((k_CD*a_CD*r) + (S_0_CD))/S_faktor
S_DE = ((k_DE*a_DE*r) + (S_0_DE))/S_faktor

% Knekningsmatrise
T_A = [1 0; 0 1; 0 1; 0 0];
T_Sym = [0; 1; -1; 0];
K_1_A = transpose(T_A)*K_1_cor*T_A;
K_G_A = transpose(T_A)*K_G_cor*T_A;
K_1_Sym = transpose(T_Sym)*K_1_cor*T_Sym;
K_G_Sym = transpose(T_Sym)*K_G_cor*T_Sym;

%SYMMETRISK KNEKNING
Lambda_Sym = -Factor_K_G / Factor_K_1;
KG_Sym = double(K_G_Sym/L^2)
K1_Sym = double(K_1_Sym/L^2)

[eigenVectors_Sym, eigenValues_Sym] = eig(K1_Sym,KG_Sym)
egenverdier_Sym = max(eigenValues_Sym)
Minste_egenverdi_Sym = min(egenverdier_Sym)
digits(5);
P_euler_sym = vpa(Minste_egenverdi_Sym / (pi^2*Lambda_Sym))
P_kr_sym = vpa(Minste_egenverdi_Sym / Lambda_Sym)


disp('husk å sette inn disse faktorene foran matrisen')
% ANTISYMMETRISK KNEKNING
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
elseif length(i_3) == 2
    i_3 = i_3(1)*i_3(2)
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
elseif length(j_3) == 2
    j_3 = j_3(1)*j_3(2)
end

K_1_A_factorized = K_1_A / i_3
K_G_A_factorized = K_G_A / j_3

Lambda = -Factor_K_G*j_3 / (Factor_K_1*i_3)
KG = double(K_G_A_factorized.*[1 1/L; 1/L 1/L^2])
K1 = double(K_1_A_factorized.*[1 1/L; 1/L 1/L^2])

[eigenVectors_A, eigenValues_A] = eig(K1,KG)
egenverdier_A = max(eigenValues_A)
Minste_egenverdi_A = min(egenverdier_A)
digits(5);
P_euler_A = vpa(Minste_egenverdi_A / (pi^2*Lambda))
P_kr_A = vpa(Minste_egenverdi_A / Lambda)
if abs(eigenValues_A(1)) < abs(eigenValues_A(4))
    eigenVectors_A(1)/eigenVectors_A(2)
else 
    eigenVectors_A(4)/eigenVectors_A(3)
end





%P = -100 ;
%L = 4;
%p = -24;
%sbc = vpa(subs(S_BC))
%sab = vpa(subs(S_AB))
%sbd = vpa(subs(S_CD))

% 1Dof bjelke
%S_AB_0 = p_AB*L_AB^2/8*[0; -1] + P_AB*L_AB/16*[0;-3]
% 2Dof bjelke
%R_0_AB = transpose(a_AB)*S_AB_0
%R_0_BC = transpose(a_BC)*(p_BC*L_BC^2/12*[1; -1] + P_BC*L_BC/8*[1;-1])
%R_0_CD = transpose(a_CD)*(p_CD*L_CD^2/12*[1; -1] + P_CD*L_CD/8*[1;-1])
%R_0 = (R_0_AB + R_0_BC + R_0_CD)/ R_0_faktor
%R = -R_0
% oppgave c)
%S_faktor = P*L/32
%r = P*L^2/(48*EI)*[2; -3]
%S_AB = ((k_AB*a_AB*r) + p_AB*L_AB^2/12*[1; -1] + P_AB*L_AB/16*[0;-3])/S_faktor
%S_BC = ((k_BC*a_BC*r) + p_BC*L_BC^2/12*[1; -1] + P_BC*L_BC/8*[1;-1])/S_faktor
%S_CD = ((k_CD*a_CD*r) + p_CD*L_CD^2/12*[1; -1] + P_CD*L_CD/8*[1;-1])/S_faktor


