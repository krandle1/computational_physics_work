#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 22:02:45 2019
@author: kirstenrandle
"""

import numpy as np
import matplotlib.pyplot as plt

#################################
#set paramters of potential
V0 = 1.
a = 1.
rmax = (2.**(1./6.))*a
#################################

#################################
#create a particle subject to a leonard-jones potential

class particle(object):
    def __init__(self, mass=1., x=0., y=0., vx=0., vy=0.):
        self.mass = mass
        self.x = x
        self.y = y
        self.vx = vx
        self.vy = vy

    def get_force(self): #force from the leonard jones potential
        r = np.sqrt(self.x**2+self.y**2)
        F = -4*V0*(6*a**6/r**7-12*a**12/r**13)
        fx = F*self.x/r
        fy = F*self.y/r
        return(fx,fy)

    def verlet(self,dt): #verlet approximation
        (fx,fy) = self.get_force()
        self.x += self.vx*dt + 0.5*fx*dt**2
        self.y += self.vy*dt + 0.5*fy*dt**2
        self.vx += 0.5*fx*dt
        self.vy += 0.5*fy*dt
        (fx,fy) = self.get_force()
        self.vx += 0.5*fx*dt
        self.vy += 0.5*fy*dt
###############################


#initialize particle parameters
m = 1.

Energy  = [i*V0 for i in np.logspace(-1.,2.,16)] #probe various energies

db = rmax/100 #impact paramter range
b = [i*db for i in range(300)]

dt = 2**-9 #time step size

#########################################
#prep plots
plt.figure(0)
plt.xlabel("x")
plt.ylabel("y")
plt.title(r"Trajectories of particles with varying b at $E = V_0$")
plt.ylim(-6,6)
plt.xlim(-6,6)
bcbar = plt.colorbar(plt.cm.ScalarMappable(cmap="cool"))
bcbar.set_ticks([0,1])
bcbar.set_ticklabels(['0',r'$3r_{min}$'])
bcbar.ax.set_ylabel("Initial b",rotation=270)

plt.figure(1)
plt.xlabel("b")
plt.ylabel(r"$\Theta$")
plt.axvline(rmax,color='r',ls='--')
plt.text(.9,3.5,r'$r_{max}$',rotation=90,c='r',fontsize=14)
plt.title(r"Scattering angle dependence on impact parameter at varying energy")
ecbar = plt.colorbar(plt.cm.ScalarMappable(cmap="winter"))
ecbar.set_ticks([0,1/3,2/3,1])
ecbar.set_ticklabels([r'0.1 $V_0$',r'$V_0$',r' 10 $V_0$',r'100 $V_0$'])
ecbar.ax.set_ylabel("Energy",rotation=270)
##########################################

#empty containers and iterators
theta = [[0 for i in b] for e in Energy]

E_color = iter(plt.cm.winter(np.linspace(0,1,len(Energy))))

sings = [[],[]] #starting point for probe range for orbits


########################################
#do physics
print("Starting loops, total 16, Energies from 0.1 to 100")

for e in range(len(Energy)):
    E=Energy[e]
    print("Energy = " + str(E))
    v = np.sqrt(2.*E/m) #using defined energy, find velocity of particle
    if E == V0:
        plt.figure(0)

    b_color = iter(plt.cm.cool(np.linspace(0,1,len(b))))

    for i in range(len(b)): #follow a particle through the potential
        p = particle(m,-3.*rmax,b[i],v,0.) #create your particle
        xpos = [p.x]
        ypos = [p.y]
        t = [0]
        j=0
        while t[j]<10:
            j+=1
            t+=[t[j-1]+dt]
            p.verlet(dt)
            xpos+=[p.x]
            ypos+=[p.y]
        if E == V0:
            bc = next(b_color)
            plt.plot(xpos,ypos,c=bc)

        theta[e][i] = np.arctan(p.vy/p.vx) #compute final angle
        if p.vx < 0:
            theta[e][i] += np.pi #pay attention to quadrant of angle

    plt.figure(1)
    ec = next(E_color)
    plt.plot(b,theta[e],c=ec)

    if max(theta[e])>np.pi: #find probe range for orbits at relevant energies
        sings[0] += [e]
        sings[1] += [b[next(a for a, v in enumerate(theta[e]) if v > np.pi)]]    
##################################

print("Loops finished")
#################################
#compute differential cross section, approaching but skipping theta = 0
#to avoid dividing by zero which will cause an overflow
dbdth = [[abs(2*db/(theta[i][j+1]-theta[i][j-1])) for j in range(1,98)] for i in range(len(theta))]
bsinth =[[b[j]/np.sin(theta[i][j]) for j in range(1,98)] for i in range(len(theta))]
diffxsec = [[bsinth[i][j]*dbdth[i][j] for j in range(len(dbdth[0]))] for i in range(len(dbdth))]



plt.figure(2)
C_color = iter(plt.cm.rainbow(np.linspace(0,1,len(theta))))
for i in range(len(theta)):
    c=next(C_color)
    plt.plot(theta[i][1:98],diffxsec[i],color=c)
plt.ylabel(r'$\frac{d\sigma}{d\Omega}$',rotation=0)
plt.xlabel(r'$\Theta$')
plt.ylim(0,500)
plt.xlim(0,np.pi)
plt.title('Differential cross section for scattering angles at varying energy')
ccbar = plt.colorbar(plt.cm.ScalarMappable(cmap="rainbow"))
ccbar.set_ticks([0,1/3,2/3,1])
ccbar.set_ticklabels([r'0.1 $V_0$',r'$V_0$',r' 10 $V_0$',r'100 $V_0$'])
ccbar.ax.set_ylabel("Energy",rotation=270)

plt.show()

input("Press any key to exit")

exit()
