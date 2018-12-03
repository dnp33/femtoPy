from femtoPy.imports import *
import femtoPy.diffusion.classes as classes
import femtoPy.PL as pl

#############
'define grid'
#############
dt=0.001
dy=0.2
y_min=0
y_max=20
t_min=0
t_max=10

G=classes.grid(dt=dt,dy=dy,y_min=y_min,y_max=y_max,t_min=t_min,t_max=t_max)

########################
'define material params'
########################
N0=1.05
s=8.5 # 8.5e10
A=1/2.1
B=1./6
C=1./2.1**3
D=1
q=-1
tMax=375
tMin=375
alpha=0.22
mu=8.5
mstar=0.067

mat=classes.material(N0=N0,s=s,A=s,B=B,C=C,D=D,q=q,tMax=tMax,tMin=tMin,alpha=alpha,mu=mu,mstar=mstar)
mat2=classes.material(N0=N0,s=s,A=s,B=B,C=C,D=D,q=q,tMax=tMax,tMin=tMin,alpha=alpha,mu=mu,mstar=mstar)

########################################
'define initial distribution and calc D'
########################################
alpha0=2.0
d0=np.exp(-alpha0*G.depth())
#d0[0]=0
T=200*np.exp(-G.time()/0.3)+300
T2=250*np.exp(-G.time()/.3)+300

dist=classes.distribution(G,d0=d0,T=T)
dist2=classes.distribution(G,d0=d0,T=T2)

##################
'container for PL'
##################
E=np.linspace(1.35,1.7,100)
gamma=0.005
Eg=1.42
theta=1
alpha0=5

PL=pl.PLspec(G,E=E,gamma=gamma,Eg=Eg,alpha0=alpha0,theta=theta,dist=dist)
PL2=pl.PLspec(G,E=E,gamma=gamma,Eg=Eg,alpha0=alpha0,theta=theta,dist=dist)



#####################
'old update function'
#####################
def update_1Iter(G,dist,mat,PL):
    f1=RHS_function(G,dist,mat)
    RHS=diff.calc_RHS(G,f1,dist)
    dist.density[:,dist.i]=diff.solve(G,f1,RHS)
    diff.calc_D(e_dens,e_params)
    f2=RHS_function(G,e_dens,e_params)
    e_dens.density[:,e_dens.i]=diff.solve(G,f2,RHS)
    diff.calc_D(dist,mat)
    PL.calcPL(dist)

