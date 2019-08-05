import Multilayer as mlt
from matplotlib.pyplot import subplots, show
from numpy import linspace, empty, sqrt, absolute

n=2

l=[1,2,4]; eps_r=[4,4,4]

# print(mlt.multiLayer.__init__.__doc__)
k0=linspace(0.001,0.01,2)
t=empty(k0.size,dtype=complex)

nTrans=n
epsTrans=nTrans*nTrans

nRefl=n
epsRefl=nRefl*nRefl
for i in range(k0.size):
    S=mlt.multiLayer(l=l,eps_r=eps_r,k0=k0[i],epsTrans=epsTrans,epsRefl=epsRefl)
    t[i]=S.S12[0,0]
print(t)
fig,ax=subplots(figsize=(5,4))
ax.plot(k0,absolute(t),color='k')
ax.axhline(absolute(2*nTrans/(nRefl+nTrans)),color='C0',dashes=(8,5))
ax.axhline(absolute(2*nRefl/(nRefl+nTrans)),color='C1',dashes=(8,5))

#ax.set_ylim(0,2)
show()
