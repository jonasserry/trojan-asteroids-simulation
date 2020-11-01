'''
Created on 26 Apr 2019

@author: jonasserry
'''
from Core import *


OrbSys = OrbitalSystem(params)

fig = plt.figure()
ax = fig.add_subplot(111)


plt.plot(OrbSys.L4[0], OrbSys.L4[1], color='black', marker='+',markersize=8)
plt.text(OrbSys.L4[0]+0.2,OrbSys.L4[1],"L4")

plt.plot(OrbSys.L4[0], -OrbSys.L4[1], color='black', marker='+',markersize=8)
plt.text(OrbSys.L4[0]+0.2,-OrbSys.L4[1]-0.3,"L5")

plt.plot(-OrbSys.rs,0,"yo", markersize = 20)
plt.text(-OrbSys.rs+0.5, 0.3, "Sun")

plt.plot(OrbSys.rp,0,"ro",markersize = 10)
plt.text(OrbSys.rp+0.2,0.3,"Jupiter")



plt.ylim([-6,6])
plt.xlim([-6,6])

for direction in ["left","bottom"]:
    ax.spines[direction].set_position('zero')
    ax.spines[direction].set_smart_bounds(False)
for direction in ["right","top"]:
    ax.spines[direction].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

circle = plt.Circle((0, 0), 5.2, color="red", fill=False)
ax.add_artist(circle)

plt.xlabel("X position /Au",position = (0.7,-0.5),size=9)
plt.ylabel ("Y position /Au",position = (-0.2,0.8),rotation = "vertical",size=9)

plt.show()