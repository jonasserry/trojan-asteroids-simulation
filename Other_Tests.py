'''
Created on 14 Apr 2019

@author: jonasserry
'''
from Core import *


#Initialise
N_Orbits = 100
Res = 100 #100 Samples per orbit


params = [G,1,0.001,R]   
OrbSys = OrbitalSystem(params)
y0 = [OrbSys.L4[0],OrbSys.L4[1],0,0]
O = orbit(OrbSys,y0)
