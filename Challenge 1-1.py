# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 12:41:05 2019

@author: kirstenrandle
"""

##Auth. Prof. Adrian Feiguin
##edited by Kirsten Randle

import numpy as np
import matplotlib.pyplot as plt


#initialize properties of system
T0 = 10.   # initial temperature
Ts = 83.   # temp. of the environment
r = 0.1    # cooling rate
tmax = 10. # maximum time


#define approximation of derivative
def euler(y,f,dx):
    approx = y + f*dx
    return approx

#initialize arrays
time = []
dt = [1.]  # initial time step
labels = ["Exact soln"]

nsteps = int(tmax/dt[0])

#compute exact result for comparison
time_exact = [x*dt[0] for x in range(nsteps+1)]
temp_exact = [Ts+(T0-Ts)*np.exp(-r*t) for t in time_exact]
plt.plot(time_exact,temp_exact, label = labels[0])

#use computed temperature at 10 seconds as benchmark for convergence
T_10 = [temp_exact[-1]] 

print ("Tempurature at 10s")

j=0
print ('{:s}  {:20.18f}'.format(labels[j],T_10[j]))
c = 1

#vary time step size and compute value at 10s
while c > 0.001: #stop when improvement from previous step size is less than 0.1%

    nsteps = int(tmax/dt[j])  # number of steps
    time = [x*dt[j] for x in range(nsteps+1)] #for plotting
    temp_e = [T0] #initialize temperature

    #iteratively compute temperature for 10s
    for i in range(1,nsteps+1):
        temp_e += [euler(temp_e[-1],-r*(temp_e[-1]-Ts),dt[j])]
    T_10 += [temp_e[-1]]

    dt += [dt[j]/2.] #pick next step size

    labels += ["dt = "+str(dt[j])+ "s"]
    plt.plot(time,temp_e,label=labels[j+1])
    print ('{:s}  {:20.18f}'.format(labels[j+1],T_10[j+1]))

    c = abs((T_10[-2]-T_10[-1])/T_10[-2])*100 #compute percent improvement
    j+=1

plt.legend(loc='best')
plt.xlim((8,10))
plt.ylim((50,60))
plt.xlabel("Time(s)")
plt.ylabel("Temperature (C)")

plt.title("Accuracy of Euler approx. with varying time step")
