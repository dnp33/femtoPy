from femtoPy.preamble import *
from femtoPy import diffusion as diff
from femtoPy import light
import time as time
from multiprocessing import Pool

plt.rcParams['axes.grid']=False

def update(dist):
    dist.prep()
    dist.dif()
    dist.bi()
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

grid=diff.grid(dt=dt,dy=dy,y_min=y_min,y_max=y_max,t_min=t_min,t_max=t_max)
y=np.asarray(grid.y)[:,0]
dy=y[1]-y[0]
t=np.asarray(grid.t.T)[:,0]

'Material Parameters'
s=8.5
ue=3.88
e=1.6e-19
B=.5
A=1./100.
alpha0=1.62
Imax=14

'CALCULATE PL SPECTRUM'

'absorption coefficient and initial spectrum'
E=np.linspace(1.42,1.55,40)
Eg=1.42
alpha=light.sqrtDOS(E,Eg)
T=0.024
PL0=light.boltzmannPL(E,alpha,T)

def process(I0):
    'print step'
    Time=time.time()-t0
    minutes=str(int(Time/60))
    seconds=str(Time%60)[:4]
    name='I0='+str((I0))+'\ntime='+minutes+' m '+seconds+' s\n'
    print(name)

    'calculate initial density'
    dSat=light.saturation(I0,alpha=alpha0,Imax=Imax,d=y)[1]

    'initialize distribution and PL array'
    Temp=200*np.exp(-t/2.)+300
    Sat=diff.dist(grid=grid,d0=dSat,s=s,u=ue,A=A,B=B,q=-e,T=Temp)
    PL=np.zeros((PL0.size,t.size))

    'loop over time steps'
    for i in range(0,grid.t.size-1):
        update(Sat)
        dens=np.asarray(Sat.density[:,i])[:,0]
        PL[:,i]=light.PL_reabsorption(alpha,PL0,dens*dens,dy)
    dens=np.asarray(Sat.density[:,-1])[:,0]
    PL[:,-1]=light.PL_reabsorption(alpha,PL0,dens*dens,dy)
    
    return

'loop over all intial intensities'
fluence=np.linspace(1,2,2)
t0=time.time()
if __name__ == '__main__':
    p = Pool(1)
    p.map(process,fluence)
    
    Time=time.time()-t0
    minutes=Time/60
    seconds=Time-Time%60
    print('\ntotal time: ',minutes,' m ',seconds,' s\n')
