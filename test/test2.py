from femtoPy.preamble import *
from femtoPy import diffusion as diff

import time as time

def update(dist):
    dist.prep()
    dist.dif()
    dist.boundary()
    dist.step()

    return

'Grid Parameters'
dt=.01
dy=.05
y_min=0
y_max=10
t_min=0
t_max=10

'Material Parameters'
s=0
ue=.85
uh=.04
T=np.zeros(np.arange(t_min,t_max+dt,dt).size)+300
e=1.6e-19
A=1./2100.


def exp(x,alpha):
    y=10*np.exp(-x*alpha)
    return y

grid=diff.grid.Grid(dt=dt,dy=dy,y_min=y_min,y_max=y_max,t_min=t_min,t_max=t_max)
e1=diff.distribution.Distribution(grid=grid,d0=np.exp(-grid.y*1.62),s=s,u=ue,A=A,q=-e,T=T)

t0=time.time()
for i in range(0,grid.t.size-1):
    update(e1)
print(time.time()-t0)

t=np.asarray(grid.t.T)[:,0]
y=np.asarray(grid.y)[:,0]
TOT=np.sum(e1.density[:,0])
rho=np.asarray(np.sum(e1.density,axis=0).T)[:,0]/TOT
fig,ax=plt.subplots(figsize=(10,8))
ax.plot(t,rho)

plt.show()
