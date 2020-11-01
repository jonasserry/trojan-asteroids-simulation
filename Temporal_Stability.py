'''
Explores the evolution of Orbits for long periods under different initial
conditions and evaluates changes in stability over time
Note: Configured to run in some specifed state. 
'''

from Core import *

"""------------------------Set Up System----------------------------"""

DefaultNumOrbits = 30
DefaultResolution = 50 #Samples per orbit

params = planets["Jupiter"]

OrbSys = OrbitalSystem(params)

"""
---------------------------------------------------------------------
                        Predefined Functions
---------------------------------------------------------------------
"""

def exploreTemporalStability(Range,spatialRes,MaxOrbits,stepSize=10,loc = OrbSys.L4, Evaluator = MaxAngularLibration):
    """Returns wander for multiple orbits along a range of initial radial positions with given spatial resolution """
    
    Alpha = np.arctan2(loc[1],loc[0]) #Initial Angle from Origin
    r = sqrt(loc[0]**2+loc[1]**2) #Radius from origin
    
    radii = np.linspace(Range[0], Range[1],spatialRes)
    
    Orbits = np.arange(MaxOrbits,step = stepSize)
    Wanders = np.full((len(radii),len(Orbits)),None)
    
    i=0
    
    for i in range(len(radii)):
        x = radii[i]
        O = orbit(OrbSys,[x*cos(Alpha),x*sin(Alpha),0,0])
        
        for j in range(len(Orbits)): #Evolve orbits in given steps
            try:
                O.evolve(stepSize,20,"odeint")
            except AssertionError:
                print("Integrator Failed")
                break
            W = Evaluator(O,O.initialPosition) #Not  efficient as re-evaluates entire array each time but also not limiting step
            Wanders[i][j] = W
            j+=1
    
        i+=1
    return (radii,Wanders,Orbits)

def numOrbitsUntillThreshold(Range,threshold,MaxOrbits,spatialRes,stepSize=10,loc = OrbSys.L4, Evaluator = MaxAngularLibration):
    """Returns number of orbits required to reach a given stability threshold. If not reached then set to MaxOrbits"""
    Alpha = np.arctan2(loc[1],loc[0]) #Initial Angle from Origin
    r = sqrt(loc[0]**2+loc[1]**2) #Radius from origin
    
    points = np.linspace(Range[0], Range[1],spatialRes)
    
    OrbitsUntillInstability = []
    i=0
    
    while i < len(points):
        x = points[i]
        MaxW = 0 
        N = 0
        O = orbit(OrbSys,[x*cos(Alpha),x*sin(Alpha),0,0])
        while MaxW < threshold and N < MaxOrbits:
            try:
                O.evolve(stepSize,20,"odeint")
            except AssertionError:
                print("Integrator Failed")
                break
            W = Evaluator(O,O.initialPosition) #Not  efficient as re-evaluates entire array each time but also not limiting step
            MaxW = max([W,MaxW])
            N+=stepSize
        
        OrbitsUntillInstability.append(N)
        i+=1
    return (points,OrbitsUntillInstability)

        
"""
---------------------------------------------------------------------
                        -----MAIN-----                                 
---------------------------------------------------------------------
"""

"""------------------Explore Temporal Stability---------------------"""
#Run Function
radii, Wanders, Orbits = exploreTemporalStability([5.25,5.27],10,2000,stepSize=20,loc = OrbSys.L4, Evaluator = MaxAngularLibration)

#Plot
plt.figure()
for i in range(len(Wanders)):
    plt.plot(Orbits,Wanders[i]/pi, label =  "{:.4f}".format(radii[i]))
plt.xlabel("Number of Orbits")
plt.ylabel("Maximum Angular Libration / (Radians/pi)")
plt.legend(title = "Radius /AU", loc = "best")


"""---------------Explore Orbits Untill Threshold-------------------"""
RL4 = convertToPolar(OrbSys.L4[0], OrbSys.L4[1])[0]

#Run Function
Radii, NumOrbs = numOrbitsUntillThreshold([5.254,5.265],threshold = pi, MaxOrbits = 10000, spatialRes= 50, stepSize=10,loc = OrbSys.L4, Evaluator = MaxAngularLibration)

#Plot
plt.plot(Radii,NumOrbs, '-x', color = "black")
plt.xlabel("Radius / AU")
plt.ylabel("# Orbits Until Threshold Reached")


plt.show()