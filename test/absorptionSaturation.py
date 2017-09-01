from femtoPy.preamble import *
from femtoPy import diffusion as diff
from femtoPy.absorption import saturationModel as sat

plt.rcParams['axes.grid']=False

import time as time


def update(dist):
    dist.prep()
    dist.dif()
    dist.mono()
    dist.boundary()
    dist.step()

    return

'Grid Parameters'
dt=.01
dy=.02
y_min=0
y_max=20
t_min=0
t_max=10

'Material Parameters'
s=.2
ue=.85
uh=.04
T=np.zeros(np.arange(t_min,t_max+dt,dt).size)+300
e=1.6e-19
B=1./2100.
A=1./2.1

grid=diff.grid.Grid(dt=dt,dy=dy,y_min=y_min,y_max=y_max,t_min=t_min,t_max=t_max)

alpha=1.62
Imax=5

Max=0

N=10
cmap=iter(plt.cm.rainbow(np.linspace(0,1,N)))
for I0 in np.linspace(1,10,N,dtype=float):
    c=next(cmap)
    y=np.asarray(grid.y)[:,0]
    t=np.asarray(grid.t.T)[:,0]
    dSat=sat.absorption(I0,alpha=alpha,Imax=Imax,d=y)
    dExp=I0*alpha*np.exp(-y*alpha)

    Sat=diff.distribution.Distribution(grid=grid,d0=dSat,s=s,u=ue,A=A,B=B,q=-e,T=T)
    Exp=diff.distribution.Distribution(grid=grid,d0=dExp,s=s,u=ue,A=A,B=B,q=-e,T=T)
    for i in range(0,grid.t.size-1):
        update(Sat)
        update(Exp)

    np.savetxt('sat'+str(I0)+'.dat',Sat.density)
    np.savetxt('exp'+str(I0)+'.dat',Exp.density)
np.savetxt('y.dat',y)
np.savetxt('t.dat',t)
