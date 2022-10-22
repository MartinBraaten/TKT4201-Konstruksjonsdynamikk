'''
Created on 22. jan. 2019

@author: Martinskole
'''
from numpy import int

D= 0.5
f=0.05
L=392
g=9.8
dt=1

z_0 = -6
v_0 = 0
T=3

z_n= z_0
v_n=v_0
N=int(T/dt)
t=0

for n in range(N):
    z_np=z_n+dt*v_n
    v_np=v_n+dt(-(2*g*z_n/L + f*v_n*v_n/(2*D)))
    
    t +=dt #= t + dt
    print"z({0}) = {1}".format(t,z_np)
    
    z_n=z_np
    v_n=v_np