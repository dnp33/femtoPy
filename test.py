from femtoPy.imports import *
from initialize import *
import femtoPy.diffusion.functions as diff

'define update function'
def update(G,dist,mat,PL):
    f1=diff.D_d2_dx2(G,mat)-diff.B(G,mat,dist)
    #f1=diff.constD_d2_dx2(G,mat)
    #f1=diff.constBoundary(G,mat,f1)
    f1=diff.boundary(G,mat,f1)
    RHS=diff.calc_RHS(G,f1,dist)
    dist.density[:,dist.i]=diff.solve(G,f1,RHS)
    
    diff.calc_D(dist,mat)
    PL.calcPL(dist)

diff.calc_D(dist,mat)
diff.calc_D(dist2,mat2)

'run simulation'
MAX_t=G.time().size-1


        
for i in range(MAX_t):
    update(G,dist,mat,PL)
    update(G,dist2,mat2,PL2)


PL.sum()
PL2.sum()

fig,ax=figure()
ax.plot(PL.E,PL.spec)
ax.plot(PL.E,PL2.spec)

fig,ax=figure()
ax.plot(PL.E,PL2.spec-PL.spec)
ax.set_xlabel('Energy (eV)')
ax.set_ylabel(r'$\Delta$PL (a.u)')

ax.axhline(0,linestyle='--',color='k',linewidth=0.5)
ax.axvline(1.487,color='k',linewidth=0.5)

plt.show()
