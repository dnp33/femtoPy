from femtoPy.preamble import *
from femtoPy import diffusion as diff

import time as time


def update(dist  ):
    dist.prep()
    dist.dif()
    dist.boundary()
    dist.bi()
    #dist.mono()
    dist.step()

    return

'Grid Parameters'
dt=.01
dy=.05
y_min=0
y_max=20
t_min=0
t_max=10

'Material Parameters'
s=0.
ue=.85
uh=.04
T=np.zeros(np.arange(t_min,t_max+dt,dt).size)+300
e=1.6e-19
A=1./2.1
B=1.


grid=diff.grid.Grid(dt=dt,dy=dy,y_min=y_min,y_max=y_max,t_min=t_min,t_max=t_max)
e=diff.distribution.Distribution(grid=grid,d0=np.exp(-grid.y*1.62),s=s,u=ue,A=A,q=-e,T=T)

t0=time.time()
for i in range(0,grid.t.size-1):
    update(e)

fig,ax=plt.subplots(figsize=(10,8))
ax.plot(grid.y,e.density[:,-1])
print(grid.t.size)

plt.show()
