'''
Created on 31. jan. 2020

@author: Martinskole
'''


import matplotlib.pylab as plt
import numpy as np

ksi = 0.02
T_1 = 2
T_2 = 3

term_1 = 7/2 + 4/np.pi**2*np.cos(4*np.pi)
term_2 = 7/2 + 4/np.pi**2*np.cos(4*np.pi) + 4/(9*np.pi**2)*np.cos(12*np.pi)
term_3 = 7/2 + 4/np.pi**2*np.cos(4*np.pi) + 4/(9*np.pi**2)*np.cos(12*np.pi) + 4/(25*np.pi**2)*np.cos(20*np.pi)

beta1 = T_1
beta2 = 3*T_1
beta3 = 5*T_1
beta4 = T_2
beta5 = 3*T_2
beta6 = 5*T_2

teta1 = np.arctan(2*ksi*beta1/(1-beta1**2))
teta2 = np.arctan(2*ksi*beta2/(1-beta2**2))
teta3 = np.arctan(2*ksi*beta3/(1-beta3**2))
teta4 = np.arctan(2*ksi*beta4/(1-beta4**2))
teta5 = np.arctan(2*ksi*beta5/(1-beta5**2))
teta6 = np.arctan(2*ksi*beta6/(1-beta6**2))

t = np.linspace(0,2,1001)

x1 = 1/(ksi*(np.sqrt((1-beta1**2)**2+(2*ksi*beta1)**2))) 
x2 = 1/(ksi*(np.sqrt((1-beta2**2)**2+(2*ksi*beta2)**2))) 
x3 = 1/(ksi*(np.sqrt((1-beta3**2)**2+(2*ksi*beta3)**2))) 
x4 = 1/(ksi*(np.sqrt((1-beta4**2)**2+(2*ksi*beta4)**2))) 
x5 = 1/(ksi*(np.sqrt((1-beta5**2)**2+(2*ksi*beta5)**2))) 
x6 = 1/(ksi*(np.sqrt((1-beta6**2)**2+(2*ksi*beta6)**2))) 

u = 7/2 + term_1*x1*np.cos(2*np.pi*t-teta1)

plt.plot(t,u)
plt.show()




