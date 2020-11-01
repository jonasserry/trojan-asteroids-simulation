'''
Embeds core functionalities of the simulation in classes and functions. 
'''

import math
import numpy as np
from numpy import pi , sin , cos , sqrt 
import matplotlib.pyplot as plt
import scipy.integrate as spInt
import timeit


"""---------------------------------------------------------------"""
"""Essential values defined here for ease of use later"""

#Defined Constants

G = 4*pi**2 #Gravitational constant
Ms = 1  #Solar Mass
Mp = 0.001#Planetary Mass
R = 5.2  #Average Distance from sun to planet /AU

params = [G,Ms,Mp,R]

#Defines various planets and constants used
planets = {
    "Jupiter": [G,Ms,Mp,R],
    "Mars": [G,1 , 3.213*10**-7,1.524],
    "Earth": [G,1, 0.000003003, 1.00],
    }

"""
---------------------------------------------------------------------
                    Core Simulator Functionality
---------------------------------------------------------------------
"""

class OrbitalSystem(object):
    """ Class defining an orbital system within which orbits can occur
    Parameters
    ----------
    params: Basic Constants of the system (G,Ms,Mp,R)
        G: Gravitational Constant
        Ms: Mass of the larger central body (Sun)
        Mp: Mass of the smaller planet (Jupiter)
        R: Distance (in Au) of Planet and Sun
    
    """

    def __init__(self, params):
        
        self.G = params[0]
        self.Ms = params[1]
        self.Mp = params[2]
        self.R = params[3]
        self.Mu = self.Mp/(self.Ms+self.Mp)
        
        self.rs = self.Mp*self.R/(self.Mp+self.Ms) #Radius from Barycentre to Sun
        self.rp = self.Ms*self.R/(self.Mp+self.Ms) #Radius from Barycentre to Planet
        
        self.w = sqrt(4*pi**2*(self.Mp+self.Ms)/self.R**3) #Angular Velocity of rotating frame
        self.period = 2*pi/self.w
        
        
        self.L4 = [self.rp - self.R/2,self.R*sin(pi/3)] #Position of L4
        self.L5 = [self.rp - self.R/2, - self.R*sin(pi/3)]
        self.Orbits = []
        
    def f (self,t,z,args = None):  
        """Returns system of 4 coupled differential equations describing equations of motion for orbits in a rotating frame of reference. Can be called by integrator class."""

        #Redefined for convenience
        rx, ry, vx, vy = z[0],z[1],z[2],z[3]
        G = self.G
        Ms = self.Ms
        Mp = self.Mp
        rp = self.rp
        rs = self.rs
        w = self.w
        
        return [vx,
           vy,
           -G * (Ms*(rx+rs)/(((rs+rx)**2 + ry**2)**(3/2)) + Mp*(rx-rp)/(((rp-rx)**2 + ry**2)**(3/2))) + 2*w*vy + rx*w**2,
           -G * (Ms*ry/(((rs+rx)**2 + ry**2)**(3/2)) + Mp*ry/(((rp-rx)**2 + ry**2)**(3/2))) - 2*w*vx + ry*w**2]
    
class orbit(OrbitalSystem):
    """Class encapsulating a single orbit within a given Orbital System
    
     Parameters
    ----------
    System: OrbitalSystem Class defining system within which orbits take place
    y0: initial position and velocity of orbit at t=0 [x0,y0,Vx0,Vy0]

    
    Attributes:
        x: Coordinates (spatial and velocity) of computed orbital trajectory
        initialPosition: starting position and velocity of Orbit
        L4,L5: Calculated of 4th and 5th Lagrange points for specified system
        
    Functionalities:
        evolve(N,S,Integrator): evolve system for N orbits with S samples per orbit using  given integrator
        
        All 'get' functions self explanatory
    """
    
    def __init__(self,System,y0):
        self.System = System
        
        self.currentPosition = y0
        
        self.initialPosition = y0[0:2]
        
        self.x = []
        self.timesteps =[]
        self.timetaken = 0
        
        self.L4 = System.L4 #redefining for convenience
        
        System.Orbits.append(self)

    def evolve (self, N, S, Integrator = "dopri5", SolveIVP = False):
        """evolve orbit for N orbits with S samples per orbit using  given integrator"""
        
        Start_Time = timeit.default_timer()
        
        T = N*2*pi/self.System.w
        Res = S*N
        
        #Use the appropriate integrator
        if SolveIVP == True:
            sol = spInt.solve_ivp(self.System.f,(0,T),self.currentPosition,method = Integrator,vectorized = True)
            y =np.transpose(sol.y)
            self.x.append(y)
            self.currentPosition = y[len(y)-1]
            self.timesteps.extend(sol.t)

        elif Integrator == "odeint":
                TSpace = np.linspace(0.001,T,Res)
                result = spInt.odeint(self.System.f, self.currentPosition, TSpace,tfirst=True)
                self.x.append(result)
                self.currentPosition = result[len(result)-1]
                self.timesteps.extend(TSpace)
        
        else:
            I = spInt.ode(self.System.f).set_integrator(Integrator)
            I.set_initial_value(self.currentPosition)
            I.set_f_params(4)
            
            TSpace = np.linspace(0.00001,T,Res)#starting from zero results in odd errors 
            for t in TSpace:
                try:
                    assert I.successful() 
                except AssertionError:
                    raise (AssertionError)
                I.integrate(t)
                self.timesteps.append(t)
                self.x.append(I.y)
            self.currentPosition = I.y
            
        timetaken = timeit.default_timer() - Start_Time
        self.timetaken = self.timetaken + timetaken
        return ()
    
    def getPositions(self):
        x = np.transpose(np.asarray(self.x))
        return(x[0],x[1])
    
    def getVelocities(self):
        return (self.x[:][2],self.x[:][3])
    
    def getTimes(self):
        return (self.timesteps)
    
    def getTimeTaken(self):
        return(self.timetaken)
    
    def Print(self):
        print("Orbit with initial position: ", self.initialPosition)
        print("Current Position: ", self.currentPosition)

"""
---------------------------------------------------------------------
                    Some Core Functions
---------------------------------------------------------------------
"""

def convertToPolar(x,y):
    """Converts x,y positions to polar coordinates r, theta"""
    theta = np.arctan2(y,x)
    theta =(theta + 2 * np.pi) % (2 * np.pi) #Convert to 0 to 2pi() values
    r = np.sqrt(x**2 + y**2)
    return(r,theta)
    
def RadialWander(Orbit,fromLoc):
    """Returns displacement from some location for each position in an Orbit"""
    x,y = Orbit.getPositions()
    x0= fromLoc[0]
    y0= fromLoc[1]
    return np.sqrt((x-x0)**2+(y-y0)**2)

def MaxRadialWander(Orbit, fromLoc):
    """Returns max displacement from some location over the entire orbit"""
    return(np.max(RadialWander(Orbit, fromLoc)))

def AvgRadialWander(Orbit, fromLoc):
    """Returns average displacement from some location"""
    return(np.average(RadialWander(Orbit, fromLoc)))

def MaxAngularLibration(Orbit,From=None):
    """Returns maximum angular libration of a given Orbit"""
    x,y = Orbit.getPositions()
    r,theta = convertToPolar(x,y)
    return( np.max(theta) - np.min(theta))

def MaxLibration(Orbit, fromLoc):
    """Returns maximum arclength of a given orbit"""
    theta = MaxAngularLibration(Orbit, fromLoc)
    x0 = fromLoc[0]
    y0 = fromLoc[1]
    return(theta*sqrt(x0**2 + y0**2))

"""
---------------------------------------------------------------------
                    Coordinate Plotting Function
---------------------------------------------------------------------
"""

def plotSystem(OrbSys, planet = "Jupiter"):
    """Plots standard coordinate system for given OrbSys and planet"""
    plt.figure()
    
    offset = OrbSys.R / 20
    
    plt.plot(OrbSys.L4[0], OrbSys.L4[1], color='black', marker='+',markersize=8)
    plt.text(OrbSys.L4[0]+offset,OrbSys.L4[1],"L4")

    plt.plot(-OrbSys.rs,0,"yo", markersize = 20)
    plt.text(-OrbSys.rs+offset, offset, "Sun")

    plt.plot(OrbSys.rp,0,"ro",markersize = 8)
    plt.text(OrbSys.rp-offset,offset,planet)
    
    plt.xlabel("X position /Au")
    plt.ylabel ("Y position /Au")
    
    return()

