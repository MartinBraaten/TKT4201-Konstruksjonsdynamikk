'''
Created on 17. feb. 2020

@author: Martinskole
'''


import matplotlib.pylab as plt
import numpy as np
from scipy.misc.common import derivative

#1a)

ksi = 0.05
rho = 2403
E_c = 24821




#1c)
t = np.linspace(0,10,1001)

G = complex(-0.0340,0.0356)
i = complex(0,1)
y = G * np.exp(i*2*t)
u_top = y * (1-np.cos(np.pi/2))
u_mid = y* (1-np.cos(np.pi/4))

#plt.plot(t,u_top)
#plt.plot(t,u_mid)
#plt.show()

#2d)

h = np.linspace(0.2,1,1001)
beta = 0.3275*np.sqrt((1+3*h)/h**3)
M = 1
W = h**2/15

sigma_st = 225000*(1+3*h)/h**2 *10**-6

sigma_0 = 0.16/(h**2*np.sqrt((1-beta**2)**2+(2*ksi*beta)**2)) 
sigma_dyn = sigma_st + sigma_0

plt.plot(h,sigma_0)
plt.plot(h,sigma_st)
plt.plot(h,sigma_dyn)
plt.xlabel("h")
plt.ylabel("Sigma [MPa]")
plt.legend(["Dynamic edge stress", "Static edge stress","Total"], loc='best', frameon=False)

plt.show()


#2e)
L = 6
x = np.linspace(0,5,1001)
phi = np.sin(x*np.pi/L)
phi1 = -np.pi**2/L**2*phi
E = 2*10**10 *10**-6
I=1
w = 40*np.pi
M_1 = E*I*phi1*y
sigm = M_1/W * 10660*15/h**2 * 1/np.sqrt((1-beta**2)**2+(2*ksi*beta)**2)

plt.plot(h,sigm)
plt.show()

