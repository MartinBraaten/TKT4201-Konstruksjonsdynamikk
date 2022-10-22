clear all
close all
clc

EI  =  20;
L  =  5;
m =  3;
Kspring = 3;
pointMass = 1 ;

%  The c o n t r i b u t i o n s from t h e 4  beam  e l e m e n t s  a r e g i v e n   a s :

K1  =  EI/L^3*[12 -6*L 0 0 0; -6*L 4*L^2 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0] ;
K2  =  EI/L^3*[0 0 0 0 0 ; 0 4*L^2 2*L^2 0 0; 0 2*L^2 4*L^2 0 0; 0 0 0 0 0; 0 0 0 0 0] ;
K3  =  EI/L^3*[12 0 -6*L 0 0; 0 0 0 0 0; -6*L 0 4*L^2 0 0; 0 0 0 0 0; 0 0 0 0 0] ;
K4  =  EI/L^3*[12 0 6*L -12 6*L; 0 0 0 0 0; 6*L 0 4*L^2 -6*L 2*L^2; -12 0 -6*L 12 -6*L; 6*L 0 2*L^2 -6*L 4*L^2];

M1  = m*L/420*[156 -22*L 0 0 0;-22*L 4*L^2 0 0 0 ;0 0 0 0 0 ;0 0 0 0 0 ;0 0 0 0 0] ;
M2  = m*L/420*[0 0 0 0 0;0 4*L^2 -3*L^2 0 0; 0 -3*L^2 4*L^2 0 0; 0 0 0 0 0; 0 0 0 0 0] ;
M3  = m*L/420*[156 0 -22*L 0   0 ;0 0 0 0 0 ;-22*L 0 4*L^2 0 0; 0 0 0 0 0; 0 0 0 0 0] ;
M4  = m*L/420*[156 0 22*L 54 -13*L ; 0 0 0 0 0; 22*L 0 4*L^2 13*L -3*L^2; 54 0 13*L 156 -22*L; -13*L 0 -3*L^2 -22*L 4*L^2];


%  F i n i s h   t h e   s c r i p t  ,   s o   t h a t   you   o b t a i n   t h e   s y s t e m   mass?
%  and   s t i f f n e s s   m a t r i c e s  M  and  K!
%  SOLUTION
%  Need   t o   add   t h e   c o n t r i b u t i o n s   from   t h e   s p r i n g   K1 ,
%  t h e   p o i n t   mass  M1 ,   a s   w e l l   a s   t h e   mass   o f   t h e   h o r i z o n t a l
%  beam   c o n t r i b u t i n g   t o   t h e   t r a n s l a t i o n a l  DOF  r1 :

K5  =  [Kspring 0 0 0 0;0 0 0 0 0; 0 0 0 0 0 ;0 0 0 0 0; 0 0 0 0 0] ;
M5  =  [pointMass+m*L  0 0 0 0;0 0 0 0 0; 0 0 0 0 0;0 0 0 0 0; 0 0 0 0 0] ;

K  =  K1  +  K2  +  K3  +  K4  +  K5 
M =  M1  +  M2  +  M3  +  M4  +  M5 