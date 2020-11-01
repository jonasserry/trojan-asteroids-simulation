'''
Explores effect of varying mass ratios on wander
'''
from Core import *

"""-------------------------Set Up System--------------------------"""

N_Orbits = 500
Res = 50#Samples per orbit

"""
---------------------------------------------------------------------
                        Predefined Functions
---------------------------------------------------------------------
"""

def ExploreMassRatio(Min,Max,N, Evaluator):
    """Returns wander for N systems with given range of mass ratios evolved for N_Orbits"""
    timeit.default_timer()
    
    #Set Up
    Ratios = []
    Wanders = []
    MSpace = np.linspace(Min,Max,N)
    Orbits = []
    
    #Itterate through mass ratios
    for Mp in MSpace:
        OrbSys = OrbitalSystem([G,1-Mp,Mp,R])
        ratio = OrbSys.Mu
        print(ratio)
        O = orbit(OrbSys,[OrbSys.L4[0],OrbSys.L4[1],0,0]) #set to modified Lagrange Point 
        try:
            O.evolve(N_Orbits,Res,"odeint",SolveIVP = False)
        except AssertionError:
            print("Integrator Failed at: ", Mp)
            continue
         
        Orbits.append(O)
        Wanders.append(Evaluator(O,O.L4))  
        Ratios.append(ratio)
    
    print ("Process took %s seconds to run. " % (timeit.default_timer()))
    return(Ratios,Wanders,Orbits)

"""
---------------------------------------------------------------------
                            ~~~~~MAIN~~~~~                                  
---------------------------------------------------------------------
"""

#Run Function in given range
Ratios, Wanders, Orbits = ExploreMassRatio(0.00001, 0.02, 400, MaxRadialWander)

#Plot Wanders
plt.figure()
plt.plot(Ratios,Wanders, linewidth = 0.8, color = "black")  
plt.xlabel("Mass Ratio")
plt.ylabel("Libration / Radians")


plt.show()
        