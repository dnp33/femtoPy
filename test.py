from femtoPy.imports import *
from initialize import *
import femtoPy.diffusion.functions as diff

'define differential equation'
def RHS_function(G,dist,mat):
    f=diff.D_d2_dx2(G,mat)-diff.B(G,mat,dist)
    f=diff.boundary(G,mat,f)

    return f

'define update function'
def update(G,dist,mat,PL):
    f1=RHS_function(G,dist,mat)
    RHS=diff.calc_RHS(G,f1,dist)
    dist.density[:,dist.i]=diff.solve(G,f1,RHS)
    diff.calc_D(dist,mat)
    PL.calcPL(dist)

diff.calc_D(dist,mat)
print(mat.D)

'run simulation'
MAX_t=G.time().size-1
MAX_t=1
for i in range(MAX_t):
    update(G,dist,mat,PL)

PL.sum()
fig,ax=figure()
ax.plot(PL.E,PL.spec)

plt.show()
