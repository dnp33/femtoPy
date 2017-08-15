from femtoPy.preamble import *
from femtoPy import diffusion as diff

import time as time


def update(dist,E):
    dist.prep()
    dist.dif()
    dist.boundary()
    dist.eField(E)
    dist.step()

    return

'Grid Parameters'
dt=.01
dy=.05
y_min=0
y_max=20
t_min=0
t_max=40

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
h1=diff.distribution.Distribution(grid=grid,d0=np.exp(-grid.y*1.62),s=s,u=uh,A=A,q=e,T=T)
E1=diff.poisson.Field(grid=grid,E0=np.asmatrix(np.zeros(grid.y.size)).T)

t0=time.time()
for i in range(0,grid.t.size-1):
    E1.solveGauss(e1,h1)
    E1.field[:,i]=E1.field[:,i]*1e-2
    update(e1,E1)
    update(h1,E1)

print(time.time()-t0)
        
fig,ax=plt.subplots(figsize=(10,8))
ax.set_title('No Electric Field')
ax.set_xlabel('y (um)')
ax.set_ylabel('n (a.u)')

N=np.arange(0,grid.t.size,int(grid.t.size/5))
color=iter(plt.cm.rainbow(np.linspace(0,1,N.size)))
for i in N:
    c=next(color)
    ax.plot(e1.grid.y,e1.density[:,i],label=str(grid.t[0,i])+' ps',color=c,linewidth=1)
    ax.plot(h1.grid.y,h1.density[:,i],linestyle='--',color=c,linewidth=1)
ax.legend()

print('times')
print('step:     ',e1.tstep)
print('boundary: ',e1.tboundary)
print('efield:   ',e1.tefield)
print('dif:      ',e1.tdif)
print('prep:     ',e1.tprep)
print('test:',e1.test)
plt.show()

# def update(dist,E):
#     dist.prep()
#     dist.dif()
#     dist.boundary()
#     dist.eField(E)
#     dist.step()

#     return

# e2=diff.distribution.Distribution(grid=grid,d0=np.exp(-grid.y*1.62),s=s,u=ue,A=A,q=-e,T=T)
# h2=diff.distribution.Distribution(grid=grid,d0=np.exp(-grid.y*1.62),s=s,u=uh,A=A,q=e,T=T)
# E2=diff.poisson.Field(grid=grid,E0=np.asmatrix(np.zeros(grid.y.size)).T)

# t0=time.time()
# for i in range(0,grid.t.size-1):
#     E2.solveGauss(e2,h2)
#     update(e2,E2)
#     update(h2,E2)
# print(time.time()-t0)


# fig,ax=plt.subplots(figsize=(10,8))
# ax.set_title('With Electric Field')
# ax.set_xlabel('y (um)')
# ax.set_ylabel('n (a.u)')

# color=iter(plt.cm.rainbow(np.linspace(0,1,N.size)))
# for i in N:
#     c=next(color)
#     ax.plot(e2.grid.y,e2.density[:,i],label=str(grid.t[0,i])+' ps',color=c)
#     ax.plot(h2.grid.y,h2.density[:,i],linestyle='--',color=c)

# ax.legend()
# plt.show()
