from femtoPy.imports import *
import femtoPy.diffusion.dist as dist
import femtoPy.diffusion.diff as diff

def RHS_function(G,e_dens,e_params):
    f=diff.D_d2_dx2(G,e_params)-diff.B(G,e_params,e_dens)
    f=diff.boundary(G,e_params,f)
    return f

'Grid parameters'
dt=.05
dy=.1
y_min=0
y_max=10
t_min=0
t_max=5

G=dist.grid(dt=dt,dy=dy,y_min=y_min,y_max=y_max,t_min=t_min,t_max=t_max)

'Material and relaxation parameters'
T=np.zeros(np.arange(t_min,t_max+dt,dt).size)+300
s=0.1

'initial density'
alpha0=0.99
d0=np.exp(-alpha0*G.depth())
#d0[0]=0
d0[-1]=0

'distribution'
e_dens=dist.distribution(G,d0)
e_params=dist.material(T=T,s=10)
diff.calc_D(e_dens,e_params)
#e_params.D=np.zeros(G.depth().size)+0.1

# fig,ax=figure()
# ax.plot(G.depth(),e_params.D)
# ax.set_xlabel('depth (um)')
# ax.set_ylabel('diffusion coefficient (cm2/Vs???)')

# plt.tight_layout()

# fig,ax=figure()
# ax.plot(G.depth(),e_params.mu)
# ax.set_xlabel('depth (um)')
# ax.set_ylabel('scattering rate (fs???)')

def update1(G,e_dens,e_params):
    f1=RHS_function(G,e_dens,e_params)
    RHS=diff.calc_RHS(G,f1,e_dens)
    e_dens.density[:,e_dens.i]=diff.solve(G,f1,RHS)
    f2=RHS_function(G,e_dens,e_params)
    e_dens.density[:,e_dens.i]=diff.solve(G,f1,RHS)


for i in range(G.time().size-1):
    update1(G,e_dens,e_params)

fig,ax=figure()
ax.plot(G.depth(),e_dens.density[:,0])
ax.plot(G.depth(),e_dens.density[:,100])
# ax.plot(G.depth(),e_dens.density[:,10])
# ax.plot(G.depth(),e_dens.density[:,20])
# ax.plot(G.depth(),e_dens.density[:,100])
# ax.plot(G.depth(),e_dens.density[:,200])

e_dens.i=0

def update2(G,e_dens,e_params):
    f1=RHS_function(G,e_dens,e_params)
    RHS=diff.calc_RHS(G,f1,e_dens)
    e_dens.density[:,e_dens.i]=diff.solve(G,f1,RHS)
    
for i in range(G.time().size-1):
    update2(G,e_dens,e_params)

ax.plot(G.depth(),e_dens.density[:,100],linestyle='--')
ax.set_xlim(-.2,10.2)

plt.tight_layout()

plt.show()
