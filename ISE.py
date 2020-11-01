'''
Runs and plots single orbits with given initial conditions. 
Also investigates the stabilities of multiple integrators
'''

from Core import *

"""
---------------------------------------------------------------------
                    Test Stability Single Orbit
---------------------------------------------------------------------
"""

#Initialise
N_Orbits = 100
Res = 100 #100 Samples per orbit


params = planets["Jupiter"] 
OrbSys = OrbitalSystem(params)
y0 = [OrbSys.L4[0],OrbSys.L4[1],0,0]
O = orbit(OrbSys,y0)

#Perform Orbit with specified integrator
O.evolve(N_Orbits,Res,"odeint")
x,y = O.getPositions()
t = O.getTimes()/OrbSys.period

print ("Integrator took %s seconds to run. " % (O.getTimeTaken()))

#Output Results
plotSystem(OrbSys)
plt.plot(x,y,label = "Asteroid path",linewidth = 0.5) 
plt.legend(loc = "best")

plt.figure()
plt.xlabel("# Orbits")
plt.ylabel("Radial Wander /Au")
plt.plot(t, RadialWander(O,OrbSys.L4))


"""
---------------------------------------------------------------------
                Test Stability Different Integrators
---------------------------------------------------------------------
"""

#Initialise new System
OrbSys = OrbitalSystem(params)
y0 = [OrbSys.L4[0]+0.2,OrbSys.L4[1]+0.2,0.0,0]

N_Orbits = 1000
Res = 100 #100 Samples per orbit

ListIntegrators = ["lsoda","vode","dopri5","odeint"]
ListIntegratorsIVP = ["RK45","RK23","Radau","BDF","LSODA"]

#Run Solvers
data1 = []
for Integrator in ListIntegrators:
    print(Integrator)
    O = orbit(OrbSys,y0)
    O.evolve(N_Orbits,Res,Integrator)
    data1.append([O.getTimeTaken(),MaxRadialWander(O,O.initialPosition[0:2])])
data1 = np.array(data1)

data2 = []
for Integrator in ListIntegratorsIVP:
    print(Integrator)
    O = orbit(OrbSys,y0)
    O.evolve(N_Orbits,Res,Integrator,SolveIVP = True)
    data2.append([O.getTimeTaken(),MaxRadialWander(O,O.initialPosition[0:2])])
data2 = np.array(data2)


#ODE Plot
fig, ax1 = plt.subplots()

index = np.arange(4)
bar_width = 0.35
opacity = 0.4

rects1 = ax1.bar(index, data1[:,0], bar_width,
                alpha=opacity, color='b',
                label="Time")

ax2 = ax1.twinx()
rects2 = ax2.bar(index + bar_width,data1[:,1], bar_width,
                alpha=opacity, color='r',
                label='Wander')

ax1.set_xlabel('Integrator (ODE Suite)')
ax1.set_ylabel('Time to Evolve 1000 Orbits /s')
ax2.set_ylabel('Wander /AU')
ax1.set_xticks(index + bar_width / 2)
ax1.set_xticklabels(ListIntegrators)

ax1.legend(loc = 2)
ax2.legend(loc = 1)
fig.tight_layout()


#SolveIVP Plot
fig, ax1 = plt.subplots()

index = np.arange(5)
bar_width = 0.35
opacity = 0.4

rects1 = ax1.bar(index, data2[:,0], bar_width,
                alpha=opacity, color='b',
                label="Time")

ax2 = ax1.twinx()
rects2 = ax2.bar(index + bar_width,data2[:,1], bar_width,
                alpha=opacity, color='r',
                label='Wander')

ax1.set_xlabel('Integrator (IVP Suite)')
ax1.set_ylabel('Time to Evolve 1000 Orbits /s')
ax2.set_ylabel('Wander /AU')
ax1.set_xticks(index + bar_width / 2)
ax1.set_xticklabels(ListIntegratorsIVP)

ax1.legend(loc = 2)
ax2.legend(loc = 1)
fig.tight_layout()

plt.show()

