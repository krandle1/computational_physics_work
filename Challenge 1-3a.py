# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 15:03:25 2019

@author: kirstenrandle
"""
#Code written by Prof. Adrian Feiguin
#Edited by Kirsten Randle

import numpy as np
import matplotlib.pyplot as plt
from operator import itemgetter

#define particle class to use later
class particle(object):
    
    def __init__(self, mass=1., y=0., v=0.):
        self.mass = mass
        self.y = y
        self.v = v
        
    def euler(self, f, dt):
        self.y = self.y + self.v*dt
        self.v = self.v + f/self.mass*dt
        
    def rk4(self, f, dt):
        self.v = self.v + f/self.mass*dt
        self.y = self.y + self.v*dt

#initialize physical parameters
g = 9.8            # g acceleration
mass = 0.01        # mass of the particle
y0 = [1.*10**5]          # initial position
v0 = 0.            # initial velocity
vt = 30.           # terminal velocity
k2 = g*mass/vt**2
Re = 6.37*10**6

def Grav(m,y):
    return g*m/(1-y/Re)**2
def Drag(m,v):
    return k2*v*abs(v)
    


dv = [0]
over = False

p=0

while p< 0.01:
    p1 = particle(mass, y0[-1], v0) #particle for g varying by height
    p2 = particle(mass, y0[-1], v0) #particle for constant g
    y1 = [y0[-1]] # since we do not know the size of the arrays, we define first a python list
    y2 = [y0[-1]]
    v1 = [v0] # the append method is more efficient for lists than arrays
    v2 = [v0]
    t1 = [0.]
    t2 = [0.]
        
    dt1 = [np.sqrt(2*10/g)]
    dt2 = [np.sqrt(2*10/g)]
    
    while p1.y>0.:
        p1.euler(-Grav(p1.mass,p1.y),dt1[-1])
        y1 += [p1.y]
        v1 += [p1.v]
        t1 += [t1[-1]+dt1[-1]]
        dt1 += [abs(10/p1.v)]
    while p2.y>0.:
        p2.euler(-p2.mass*g,dt2[-1])
        y2+=[p2.y]
        v2 +=[p2.v]
        t2 += [t2[-1]+dt2[-1]]
        dt2 += [abs(10/p2.v)]

    p = abs((v1[-1]-v2[-1])/v2[-1])
    y0+=[y0[-1]+100]
    
plt.plot(v1,y1,label='Variable Gravity')
plt.plot(v2,y2,label='Constant Gravity')
plt.legend(loc='best')
plt.ylabel("Height (m)")
plt.xlabel("Velocity (m/s)")
print("Initial Drop height: {:6.0f} m").format(y0[-1])