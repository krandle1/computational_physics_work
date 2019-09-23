# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 10:53:55 2019

@author: kirstenrandle
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ColorConverter as cc

class particle2(object):
    
    def __init__(self, mass=1., x=0., y=0., vx=0., vy=0.):
        self.mass = mass
        self.x = x
        self.y = y
        self.vx = vx
        self.vy = vy
       
    def euler(self, fx, fy, dt):
        self.vx = self.vx + fx/self.mass*dt
        self.vy = self.vy + fy/self.mass*dt
        self.x = self.x + self.vx*dt
        self.y = self.y + self.vy*dt
        
g = 9.8            # g acceleration
v0 = 30.           # initial velocity

dt = 0.01           # time step
c  = 0.1           # k/m    

colors = ['red','orange','yellow','green','magenta','cyan','blue','purple','black']
maxx = [0]


b=1
angle = 1
while b>=0:
    x = [0]                                  # we need to initialize the arrays for each value of the angle
    y = [0]
    vx = [np.cos(angle*0.05*np.pi/2.)*v0] 
    vy = [np.sin(angle*0.05*np.pi/2.)*v0] 
    t = [0.]
    
    p = particle2(1., 0., 0., vx[0], vy[0])
    
    while p.y >= 0.:
        fy = -p.mass*g -c*p.mass*p.vy*abs(p.vy)
        fx = -c*p.mass*p.vx*abs(p.vx)
        p.euler(fx, fy, dt)
        x.append(p.x)
        y.append(p.y)
        vx.append(p.vx)
        vy.append(p.vy)
        t.append(t[-1]+dt)
        
    t_data = np.array(t) # we convert the list into a numpy array for plotting
    x_data = np.array(x)
    y_data = np.array(y)
    vx_data = np.array(vx)
    vy_data = np.array(vy)
    maxx += [x_data[-1]]
    
   
    print("Distance at {:3.1f} Degrees: {:3.3} m").format(np.degrees(angle*0.05*np.pi/2.),x_data[-1])
    b = maxx[-1]-maxx[-2]
    angle += 1 
b=1
angle += -3
maxx=maxx[:-3]
while b>=0:
    angle += 0.5 
    x = [0]                                  # we need to initialize the arrays for each value of the angle
    y = [0]
    vx = [np.cos(angle*0.05*np.pi/2.)*v0] 
    vy = [np.sin(angle*0.05*np.pi/2.)*v0] 
    t = [0.]
    
    p = particle2(1., 0., 0., vx[0], vy[0])
    
    while p.y >= 0.:
        fy = -p.mass*g -c*p.mass*p.vy*abs(p.vy)
        fx = -c*p.mass*p.vx*abs(p.vx)
        p.euler(fx, fy, dt)
        x.append(p.x)
        y.append(p.y)
        vx.append(p.vx)
        vy.append(p.vy)
        t.append(t[-1]+dt)
        
    t_data = np.array(t) # we convert the list into a numpy array for plotting
    x_data = np.array(x)
    y_data = np.array(y)
    vx_data = np.array(vx)
    vy_data = np.array(vy)
    maxx += [x_data[-1]]
    
   
    print("Distance at {:3.1f} Degrees: {:3.3} m").format(np.degrees(angle*0.05*np.pi/2.),x_data[-1])
    b = maxx[-1]-maxx[-2]

angle+= -0.5

print("Maximum occurs at {:3.3} degrees").format(np.degrees(angle*0.05*np.pi/2.))
