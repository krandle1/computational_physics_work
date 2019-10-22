import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from IPython.display import HTML
import math
import itertools as itool

class BoundaryCondition:
    RBC, PBC = range(2)
    
class StartConf:
    Triangular, Rectangular, Random, Confined = range(4)

class particle2(object):

    def __init__(self, mass=1., x=0., y=0., vx=0., vy=0.):
        self.mass = mass
        self.x = x
        self.y = y
        self.vx = vx
        self.vy = vy
        self.fx = 0.
        self.fy = 0.
       
    def euler(self, fx, fy, dt):
        self.vx = self.vx + self.fx/self.mass*dt
        self.vy = self.vy + self.fy/self.mass*dt
        self.x = self.x + self.vx*dt
        self.y = self.y + self.vy*dt
        
    def accel(self, dt):
        self.vx = self.vx + self.fx/self.mass*dt
        self.vy = self.vy + self.fy/self.mass*dt
        
    def move(self, dt, dt2half):
        self.x = self.x + self.vx*dt + self.fx/self.mass*dt2half
        self.y = self.y + self.vy*dt + self.fy/self.mass*dt2half
class MDsystem(object):

    def __init__(self, lx, ly, N, dt, bc): 
        self.N = N
        self.lx = ly
        self.ly = lx
        self.dt = dt
        self.dt2 = dt*dt
        self.dthalf = dt * 0.5
        self.dt2half = self.dt2 * 0.5
        self.bc = bc
        self.particles = [particle2()]
        for i in range(1,N):
            self.particles.append(particle2()) # we create a list of N particles

    def init(self, config, vmax):
        np.random.seed(1272121) # arbitrary seed
        nrows = int(math.sqrt(self.N)*float(self.ly)/float(self.lx))
        ncols = int(math.sqrt(self.N)*float(self.lx)/float(self.ly))
        ax = float(lx) / float(ncols)
        ay = float(ly) / float(nrows)
        i = 0

        if config == StartConf.Triangular:
            #Challenge
            nrows = int(math.sqrt(self.N))
            ncols = int(nrows)
            ax = float(lx)/float(ncols)
            ay = float(ly)/float(nrows)
            for row, col in itool.product(range(nrows),range(ncols)):
                self.particles[i].x = col*ax+ax/2.
                if row%2 == 1:
                    self.particles[i].x += ax/2.
                self.particles[i].y = row*ay+ay/2.
                i+=1
        elif config == StartConf.Rectangular:
            for row, col in itool.product(range(nrows),range(ncols)):
                if i >= self.N: 
                    break
                self.particles[i].x = col*ax+ax/2.
                self.particles[i].y = row*ay+ay/2.
                i+=1

            for row, col in itool.product(range(1,nrows),range(1,ncols)):
                if i >= self.N: 
                    break
                self.particles[i].x = col*ax+ax/2.+ax/4.
                self.particles[i].y = row*ay+ay/2.+ay/4.
                i+=1
                
        elif config == StartConf.Confined:
            ax /= 2.
            for row, col in itool.product(range(nrows),range(ncols)):
                if i >= self.N: 
                    break
                self.particles[i].x = col*ax+ax/2.
                self.particles[i].y = row*ay+ay/2.
                i+=1

            for row, col in itool.product(range(nrows),range(ncols)):
                if i >= self.N: 
                    break
                self.particles[i].x = col*ax+ax/2.+ax/4.
                self.particles[i].y = row*ay+ay/2.+ay/4.
                i+=1
                
        elif config == StartConf.Random:
            for i in range(self.N):
                overlap = True;
                while overlap:
                    overlap = False;
                    self.particles[i].x = np.random.random()*self.lx
                    self.particles[i].y = np.random.random()*self.ly
                    for j in range(i):
                        r12 = self.distance(self.particle[i], self.particle.p[j])
                        if r12 < 1.: 
                            overlap = True;
                            break
                                
        # Velocities
        for p in self.particles:
            p.vx = vmax*(2. * np.random.random() - 1);
            p.vy = vmax*(2. * np.random.random() - 1);

        # We set total momentum to zero
        vxcm = 0.
        vycm = 0. # Velocity of the center of mass
        for p in self.particles:
            vxcm += p.vx;
            vycm += p.vy;
        
        vxcm /= self.N
        vycm /= self.N
        for p in self.particles:
            p.vx -= vxcm;
            p.vy -= vycm;
            
        self.forces()          
        
    def evolve(self):
        for p in self.particles:
            p.move(self.dt, self.dt2half)
            p.accel(self.dthalf)
            self.boundary(p)

        self.forces()

        for p in self.particles:
            p.accel(self.dthalf)

        
    def distance(self, p, other):
        (r12, dx, dy) = self.distance2(p, other)
        return math.sqrt(r12)

    def distance2(self, p, other):
        dx = other.x - p.x;
        dy = other.y - p.y;


        # nearest image convention
        if self.bc == BoundaryCondition.PBC:
            if abs(dx) > self.lx/2:
                dx -= dx*lx/abs(dx)
                
            if abs(dy) > self.ly/2:
                dy -= dy*ly/abs(dy)
        
        r12 = dx * dx + dy * dy
        if r12 == 0.:
            r12 = 0.000000000001
        return (r12, dx, dy)

    def force(self, p, other):  #Lennard-Jones
        (r12,dx,dy) = self.distance2(p, other)
        r2 = 1./r12
        r6 = r2 * r2 * r2
        f = 24.*r2*r6*(2.*r6-1.)
        fx = f*(-dx)
        fy = f*(-dy)
        return (fx, fy);
      
    def forces(self):
        # Compute the interaction forces between particles
        for p in self.particles:
            p.fx = 0.
            p.fy = 0.
    
        for i in range(self.N):
            p = self.particles[i]
            for j in range(i+1,self.N):
                other = self.particles[j]
                (fx, fy) = self.force(p, other)
                p.fx += fx
                p.fy += fy
                other.fx -= fx
                other.fy -= fy

    def boundary(self, p):
        if self.bc == BoundaryCondition.RBC:
            
            if p.x < 0 or p.x > lx:  
                p.x = lx - p.x%lx
                p.vx = -p.vx
            if p.y < 0 or p.y > ly:
                p.y = ly - p.y%ly
                p.vy = -p.vy

        elif self.bc == BoundaryCondition.PBC:

            if p.x < 0 or p.x > lx: 
                p.x = p.x%lx
            if p.y < 0 or p.y > ly:  
                p.y = p.y%ly

            
    def kinetic_energy(self): # Challenge
        ke = 0.
        for p in self.particles:
            ke += 0.5*p.mass*( math.pow(p.vx,2.) + math.pow(p.vy,2.) )
        return ke
    
    def pot_energy(self): # Challenge
        pe = 0.

        for n in range(self.N):
            p = self.particles[n]

            for i in range(n+1, self.N):
                q = self.particles[i]
                r = self.distance(p,q)
                pe += 4*( math.pow( (1/r), 12.) - math.pow( (1/r), 6.))
        return pe
    
    def total_energy(self):
        return self.kinetic_energy() + self.pot_energy()

lx = 8
ly = 8
N  = 12
dt = 0.001
v0 = 1.0

T = [1,2,4]

colors  = ['blue','yellow','red']
labels  = ['T = {}'.format(t) for t in T]

a = 10
b = 100

#Part a
for t in range(len(T)):

    Sc = MDsystem(lx, ly, N, dt, BoundaryCondition.PBC)
    Sc.init(StartConf.Confined, v0)
    Tavg = [0 for i in range(a)]
    time = [dt*b*(i+1/2) for i in range(a)]
    for i in range(a):
        tsum = 0
        for j in range(b):
            Tactual = 2*Sc.kinetic_energy()/Sc.N/3
            tsum += Tactual
            Sc.evolve()
        Tavg[i] = tsum/float(a)
        f = math.sqrt(T[t]/Tavg[i])
        for p in Sc.particles:
            p.vx *= f
            p.vy *= f
    plt.figure(0)    
    plt.plot(time,Tavg,color=colors[t],label =labels[t])

plt.figure(0)
plt.legend(loc='best')
plt.xlabel('Time')
plt.ylabel('Velocity')
plt.title('Part a')

#parts b & c
a = 200
b = 5

styles = ['-','--',':']
energies = ['KE','PE','Total']
elabels = [[ e + ' T = ' + str(t) for e in energies] for t in T]
for t in range(len(T)):

    Sf = MDsystem(lx, ly, N, dt, BoundaryCondition.PBC)
    Sf.init(StartConf.Confined, v0)
    Tavg = [0 for i in range(a)]
    time = [dt*b*(i+1/2) for i in range(a)]
    KE = [0 for i in range(a*b)]
    PE = [0 for i in range(a*b)]
    TE = [0 for i in range(a*b)]
    e_time = [x*dt for x in range(a*b)]
    for i in range(a):
        tsum = 0
        for j in range(b):
            KE[i*b+j] = Sf.kinetic_energy()
            PE[i*b+j] = Sf.pot_energy()
            TE[i*b+j] = Sf.total_energy()
            Tactual = 2*Sf.kinetic_energy()/Sf.N/3
            tsum += Tactual
            Sf.evolve()
        Tavg[i] = tsum/float(a)
        f = math.sqrt(T[t]/Tavg[i])
        for p in Sf.particles:
            p.vx *= f
            p.vy *= f
    plt.figure(1)
    plt.plot(time,Tavg,color=colors[t],label =labels[t])
    
    plt.figure(2)
    plt.plot(e_time,KE,color='magenta',linestyle = styles[t],label = elabels[t][0])
    plt.plot(e_time,PE,color='blue',linestyle = styles[t],label = elabels[t][1])
    plt.plot(e_time,TE,color='black',linestyle = styles[t],label = elabels[t][2])

plt.figure(1)
plt.legend(loc='best')
plt.xlabel('Time')
plt.ylabel('Velocity')
plt.title('Part b')

plt.figure(2)
plt.legend(loc='best')
plt.xlabel('Time')
plt.ylabel('Energy')
plt.title('Part c')
plt.show()


input('Press enter to exit')

exit()
