from femtoPy.imports import *
import femtoPy.diffusion.classes as dist
import femtoPy.diffusion.functions as diff

def RHS_function(G,e_dens,e_params):
    f=diff.D_d2_dx2(G,e_params)-diff.B(G,e_params,e_dens)
    f=diff.boundary(G,e_params,f)
    # f[1,:]=f[1,:]-np.sin(e_dens.n())**2

    return f

'Grid parameters'
dt=.1
dy=.1
y_min=0
y_max=10
t_min=0
t_max=10

G=dist.grid(dt=dt,dy=dy,y_min=y_min,y_max=y_max,t_min=t_min,t_max=t_max)

'Material and relaxation parameters'
T=np.zeros(np.arange(t_min,t_max+dt,dt).size)+300
s=1
B=.2
tMin=374

'initial density'
alpha0=0.99
d0=np.exp(-alpha0*G.depth())
#d0[0]=0
d0[-1]=0

'initialize classes'
e_dens=dist.distribution(G,d0)
e_params=dist.material(T=T,s=s,B=B,tMin=tMin)
diff.calc_D(e_dens,e_params)
#e_params.D=np.zeros(G.depth().size)+0.1

'plot and simulation params'
MAX_t=G.time().size-1
#MAX_t=5
t_test=4.9

def update1(G,e_dens,e_params):
    f1=RHS_function(G,e_dens,e_params)
    RHS=diff.calc_RHS(G,f1,e_dens)
    e_dens.density[:,e_dens.i]=diff.solve(G,f1,RHS)
    diff.calc_D(e_dens,e_params)
    f2=RHS_function(G,e_dens,e_params)
    e_dens.density[:,e_dens.i]=diff.solve(G,f2,RHS)
    diff.calc_D(e_dens,e_params)

for i in range(MAX_t):
    update1(G,e_dens,e_params)

fig,ax=figure()
ax.plot(G.depth(),e_dens.n_t(G,0))
ax.plot(G.depth(),e_dens.n_t(G,t_test))

'save final density for comparison'
d1=e_dens.n_t(G,t_test).copy()

'reinitialize code and change diffusion parameter'
e_dens.i=0
e_params.D=np.zeros(G.depth().size)+np.amax(e_params.D)

def update2(G,e_dens,e_params):
    f1=RHS_function(G,e_dens,e_params)
    RHS=diff.calc_RHS(G,f1,e_dens)
    e_dens.density[:,e_dens.i]=diff.solve(G,f1,RHS)
    f2=RHS_function(G,e_dens,e_params)
    e_dens.density[:,e_dens.i]=diff.solve(G,f2,RHS)

    
for i in range(MAX_t):
    update2(G,e_dens,e_params)

ax.plot(G.depth(),e_dens.n_t(G,t_test),linestyle='--')
ax.set_xlim(-.2,10.2)

fig,ax=figure()
ax.plot(G.depth(),e_dens.n_t(G,t_test)-d1)

ax.set_xlabel('depth (a.u)')
ax.set_ylabel('difference in density (a.u)')

plt.tight_layout()
plt.show()
