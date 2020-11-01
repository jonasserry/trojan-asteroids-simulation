'''
Explores spatial stability in One Dimension
Note: Configured to run in some specifed state. 
Some parameters may have to be changed to plot similair figures
'''

from Core import *

"""-----------------------Set Up System-------------------------"""
DefaultNumOrbits = 800
DefaultResolution = 50 #Samples per orbit

planet = "Jupiter"
params = planets[planet]
OrbSys = OrbitalSystem(params)

#For Convenience:
L4_R, L4_Theta = convertToPolar(OrbSys.L4[0], OrbSys.L4[1]) 
L5_R, L5_Theta = convertToPolar(OrbSys.L5[0], OrbSys.L5[1]) 

"""
------------------------------------------------------------------
                     Predefined Functions
------------------------------------------------------------------
"""

def exploreSpatialStabilityAlong1D(Range,loc,Samples,Evaluator,Dimension,N_Orbits = DefaultNumOrbits,Res = DefaultResolution, Integrator="odeint"):
    """Returns wander of orbits explored spatially along one dimension (radial or angular)
     in given range from given location for a specified number of samples"""

    timeit.default_timer()
    
    Alpha = np.arctan2(loc[1],loc[0]) #Initial Angle from Origin
    r = sqrt(loc[0]**2+loc[1]**2) #Radius from origin
    
    #Set Up
    Points = np.linspace(Range[0], Range[1],Samples)
    Wanders = []
    UsedPoints =[]
    Xs = []
    Ys = []
    i = 0
    
    while i<Samples: 
        #Set Point
        if Dimension == "Radial":
            X = (Points[i])*cos(Alpha)
            Y = (Points[i])*sin(Alpha)
        elif Dimension == "Angular":
            X = r*np.cos(Alpha + Points[i])
            Y = r*np.sin(Points[i] + Alpha)
        else:
            print("Unknown Dimension!")
            break
        
        #Initialise Orbit
        O = orbit(OrbSys,[X,Y,0,0])
        try:
            O.evolve(N_Orbits,Res,Integrator)
        except AssertionError:
            i+=1
            print("Integrator Failed at: ", (X,Y))
            continue
        
        #Save wander and point
        Wanders.append(Evaluator(O,O.initialPosition))
        Xs.append(X)
        Ys.append(Y)

        UsedPoints.append(Points[i])
        i+=1
        print("Runs left: ",Samples-i)

    print ("Process took %s seconds to run. " % (timeit.default_timer()))
    return(UsedPoints,Wanders,[Xs[:],Ys[:]])


def findStablePoints(y,threshold):
    """returns indices of minima below a certain threshold for given y range"""
    dy = np.gradient(y)
    indeces = [] 
    sign = dy[0]<=0 #is y negative
    i=1
    while (i<len(y)):
        if dy[i]<0:
            sign = True
        if dy[i]>0 and sign:
            if y[i]<threshold:
                indeces.append(i)
            sign = False
        i+=1
    return (indeces)
    
    
    
"""
---------------------------------------------------------------------
                            -----MAIN-----                                  
---------------------------------------------------------------------
"""

"""--------------------------Along Arc---------------------------"""

NumPoints = 500
Range = [-pi,pi]

#Run Function
Angles, Wanders, Coords = exploreSpatialStabilityAlong1D(Range, OrbSys.L4, NumPoints, MaxAngularLibration,"Angular")
stablePoints = findStablePoints(Wanders, 1)

#Plot Range explored and L4,L5 positions
plotSystem(OrbSys,planet)
plt.plot(Coords[0],Coords[1])
plt.plot(OrbSys.L5[0], OrbSys.L5[1], color='black', marker='+',markersize=8)
plt.text(OrbSys.L5[0]+0.2,OrbSys.L5[1],"L5")


#Plot stable positions found
plt.plot(np.take(Coords[0],stablePoints),np.take(Coords[1],stablePoints),"o",color ="black", markersize = 5)

#Plot Angle Vs Libration
plt.figure()
plt.plot(Angles, Wanders, linewidth = 0.6)
plt.xlabel("Angle from L4/ Radians")
plt.ylabel("Max Angular Libration/ Radians")

plt.plot(0,0.0,"+", color = "black")
plt.text(0.3,0,"L4")

plt.plot(-L4_Theta*2,0.0,"+", color = "black")
plt.text(-L4_Theta*2-0.3,0,"L5")


"""-------------------------Along Radius----------------------------"""

NumPoints = 200
Range = [L4_R-0.2,L4_R+0.2]
Radii, Wanders, Coords = exploreSpatialStabilityAlong1D(Range, [OrbSys.L4[0],OrbSys.L4[1]], NumPoints, MaxAngularLibration,"Radial")
stablePoints = findStablePoints(Wanders, 3)


#Plot Coordinate System
plotSystem(OrbSys)
plt.plot(Coords[0],Coords[1], linewidth = 0.6)

#Plot stable positions
plt.plot(np.take(Coords[0],stablePoints),np.take(Coords[1],stablePoints),"o",color ="black", markersize = 5)


#Plot Radii Vs Wanders
plt.figure()
plt.plot(Radii, Wanders, linewidth = 0.6)
plt.xlabel("Radius/ AU")
plt.ylabel("Angular Libration/ Radians")

#Plot Position of L4
plt.plot(L4_R,0.05,"+", color = "black")
plt.text(L4_R+0.0003,0,"L4")


plt.show()
    
    