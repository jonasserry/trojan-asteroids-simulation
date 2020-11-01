'''
Explores 2-Dimensional region in space and velocity and plots heat maps of stability
'''

from Core import *
import random
import matplotlib.tri as tri

"""-------------------------Set Up System--------------------------"""
DefaultNumOrbits = 30
DefaultResolution = 50 #Samples per orbit

params = planets["Jupiter"]

OrbSys = OrbitalSystem(params)

"""
---------------------------------------------------------------------
                        Predefined Functions
---------------------------------------------------------------------
"""
def exploreSpatialStabilityRandomly(lims1,lims2,SpatialRes,Evaluator,Integrator = "odeint", N_Orbits = DefaultNumOrbits,Res = DefaultResolution):
    """Explores a defined region within given limits randomly and returns wander
     computed after evolving for a given number of orbits"""
    timeit.default_timer()
    Xs = []
    Ys= []
    Wanders = []
    tts = []
    
    while SpatialRes > 0:
        #Find random point
        Xp = random.uniform(lims1[0],lims1[1])
        Yp = random.uniform(lims2[0],lims2[1])
    
        #Evolve Orbit
        O = orbit(OrbSys,[Xp,Yp,0,0])
        try:
            O.evolve(N_Orbits,Res, Integrator = "odeint")
        except AssertionError:
            print("Integrator Failed at: ", (Xp,Yp))
            continue

        #Compute and save wanders
        tts.append(O.getTimeTaken())
        Wanders.append(Evaluator(O, O.L4))
        Xs.append(Xp)
        Ys.append(Yp)
        SpatialRes -=1
        if SpatialRes%10 == 0: #Estimated time remaining
            print("Estimated seconds remaining: ", (int(np.average(tts)*SpatialRes)))

    print ("Process took %s seconds to run. " % (timeit.default_timer()))
    return(Xs,Ys,Wanders)


def exploreStabilityOnGrid(lims1,lims2,res1,res2,Evaluator,Integrator = "odeint", N_Orbits = DefaultNumOrbits,Res = DefaultResolution, polar = False, space = "spatial", OrbSys = OrbSys):
    """Explores a defined region in space or velocity within given limits on a grid 
    of defined resolution and returns wander computed after evolving for a given number of orbits
    Can be performed in polar or cartesian coordinates"""
    timeit.default_timer()
    tts = []
    
    if polar: #if lims are fed in polar coordinates then perform appropriate conversions
        
        if space == "spatial": #if space being explored is spatial
            r = np.linspace(lims1[0],lims1[1],res1)
            theta = np.linspace(lims2[0],lims2[1],res2)
            r, theta = np.meshgrid(r, theta)
            X = r*np.cos(theta)
            Y = r*np.sin(theta)
            
        if space == "velocity": #if velocity space is being explored
            #Create Polar Coordinates
            Vr = np.linspace(lims1[0],lims1[1],res1)
            Omega = theta = np.linspace(lims2[0],lims2[1],res2)
            Vr, Omega = np.meshgrid(Vr, Omega)
    
            #Convert from Polar to Cartesian for use
            r, theta = convertToPolar(OrbSys.L4[0],OrbSys.L4[1])
            X = (Vr*np.cos(theta) - r*Omega*np.sin(theta))
            Y = (Vr*np.sin(theta) - r*Omega*np.cos(theta))
            
    else: #If cartesian limits are fed
        x = np.linspace(lims1[0],lims1[1],res1)
        y = np.linspace(lims2[0],lims2[1],res2)
        X,Y = np.meshgrid(x, y)
    
    Wanders = np.full((res2,res1),None)
    
    #Itterate through initial conditions
    for i in range(res2):
        for j in range(res1):
            x0 = X[i][j]
            y0 = Y[i][j]
            
            #Initialise and Evolve Orbit
            if space == "spatial":
                O = orbit(OrbSys,[x0,y0,0,0])
            elif space == "velocity":
                O = orbit(OrbSys,[OrbSys.L4[0],OrbSys.L4[1],x0,y0])
            else:
                print("Undefined space type")
                return()
            try:
                O.evolve(N_Orbits,Res, Integrator = "odeint")
            except AssertionError:
                print("Integrator Failed at: ", (x0,y0))
                continue
            tts.append(O.getTimeTaken())
            Wanders[i][j] = Evaluator(O, O.initialPosition)
        print("Estimated seconds remaining: ", (int(np.average(tts)*(res2-i-1)*res1)))

    print ("Process took %s seconds to run. " % (timeit.default_timer()))
    
    if space == "velocity" and polar:#Return some additional values for velocity exploration
        return (X,Y,Vr,Omega,Wanders)
        
    return(X,Y,Wanders)


def AreaOfStability(W,threshold):
    """Finds fraction of a given matrix representing a grid that lies below a certain 
    threshold"""
    stablePoints = 0
    for i in range(len(W)):
        for j in range (len(W[0])):
            if W[i][j] < threshold:
                stablePoints +=1
    
    total = len(W)*len(W[0])
    return(stablePoints/total)

"""
---------------------------------------------------------------------
                          -----MAIN-----                                
---------------------------------------------------------------------
"""
"""--------------Explore 2D Spatial Stability Randomly------------"""

#Set Up System
NumOrbits = 20
Resolution = 10 #Samples per orbit
NumPoints = 5000
xlims = [-6,6]
ylims = [-6,6]


#Run Function
X,Y,W = exploreSpatialStabilityRandomly(xlims,ylims, NumPoints, MaxAngularLibration,
                                        N_Orbits = NumOrbits,Res = Resolution)

#Plot Coloured Scatter
plt.figure()
plt.scatter(X,Y,c = W, cmap = "plasma")
plt.xlabel("X position /Au")
plt.ylabel("Y position /Au")


#Plot Interpolated HeatMap in Cartesian:
triang = tri.Triangulation(X, Y)
plotSystem(OrbSys)
plt.tricontourf(triang, W, levels = 1000, cmap="RdYlBu",extend = "neither") 
cbar = plt.colorbar()
cbar.ax.set_ylabel('Maximum Libration /Radians')


#In Polar:
r, theta = convertToPolar(np.asarray(X), np.asarray(Y))
triang = tri.Triangulation(r, theta)

plt.figure()
plt.tricontourf(triang, W, 1000, cmap="RdYlBu") 
plt.xlabel("Radial Position /Au")
plt.ylabel("Angular Position / (Radians/pi)")
cbar = plt.colorbar()
cbar.ax.set_ylabel('Maximum Libration / Radians')



"""----------Explore Spatial Stability on Grid-------------"""

#Set Up System
lims1 = [5,5.4]
lims2 = [0.15,pi]
#Compute total area explored
Area = 2*(lims2[1]-lims2[0])*lims1[0]*(lims1[1]-lims1[0])

#Run Function
X,Y,W = exploreStabilityOnGrid(lims1, lims2, 100, 100, MaxAngularLibration, N_Orbits = 100, Res = 10, polar=True)

#Plot Cartesian
plt.figure()
plt.scatter(X,Y,c = W, cmap = "plasma")
plt.xlabel("X position /Au")
plt.ylabel("Y position /Au")

print(AreaOfStability(W,2)*Area)
plotSystem(OrbSys)
plt.contourf(X,Y, W, levels = 1000, cmap="RdYlBu")
cbar = plt.colorbar()
cbar.ax.set_ylabel('Maximum Libration /Radians')

#Plot Polar
r, theta = convertToPolar(X, Y)
plt.figure()
plt.scatter(r,theta,c = W, cmap = "plasma")
plt.xlabel("Radial Position /Au")
plt.ylabel("Angular Position/(Radians/pi)")

plt.figure()
plt.contourf(r,theta, W, 1000, cmap="RdYlBu") 
plt.xlabel("Radial Position /Au")
plt.ylabel("Angular Position / (Radians/pi)")
cbar = plt.colorbar()
cbar.ax.set_ylabel('Maximum Libration / Radians')


"""--------Explore Velocity Space Stability on Grid---------"""


#Limits Cartesian
lims1 = [-0.6,0.5]
lims2 = [-3.5,0.6]

X,Y,W = exploreStabilityOnGrid(lims1, lims2, 100, 100, MaxAngularLibration, N_Orbits = 30, Res = 10, space = "velocity",polar = False)
W = np.asarray(W)

#Plot Cartesian Velocities
plt.figure()
plt.scatter(X,Y,c = W, cmap = "plasma")
plt.xlabel("Vx Initial/ (Au/year)")
plt.ylabel("Vy Initial/ (Au/year)")

print(AreaOfStability(W,pi))
plt.figure()
plt.xlabel("Vx Initial/ (Au/year)")
plt.ylabel("Vy Initial/ (Au/year)")
plt.contourf(X,Y, W, levels = 1000, cmap="RdYlBu")
cbar = plt.colorbar()
cbar.ax.set_ylabel('Maximum Libration /Radians')

#Limits Polar
lims2 = [-0.3,0.3]
lims1 = [-3.0,0.5]

#Run Function In Polar Coordinates
X,Y,Vr,Omega,W = exploreStabilityOnGrid(lims1, lims2, 100,100, MaxAngularLibration, N_Orbits = 30, Res = 10, space = "velocity",polar = True)
W = np.asarray(W)

#Polar Plot
plt.figure()
plt.scatter(Vr,Omega,c = W, cmap = "plasma")
plt.xlabel("Initial Angular Velocity/ (Au/year)")
plt.ylabel("Initial Radial Velocity/ (Radians/year)")

plt.figure()
plt.contourf(Vr,Omega, W, 1000, cmap="RdYlBu") 
plt.xlabel("Initial Angular Velocity/ (Radians/year)")
plt.ylabel("Initial Radial Velocity/ (Au/year)")
cbar = plt.colorbar()
cbar.ax.set_ylabel('Maximum Libration / Radians')

plt.show()

