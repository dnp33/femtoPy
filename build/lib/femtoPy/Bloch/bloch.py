from femtoPy.imports import *

def Omega(t):
    return np.exp(-t*t)

def sBloch(P,fc,fv,det,t):
    tau1=2
    tau2=1
    dP=-1j*det-1j*Omega(t)*(fc-fv)-P/tau2
    dfc=Omega(t)*np.imag(P)-fc/tau1
    dfv=-dfc

    return dP,dfc,dfv

def rk4(P0,fc0,fv0,det,t,dt):
    dP1,dfc1,dfv1=sBloch(P0,fc0,fv0,det,t)
    dP2,dfc2,dfv2=sBloch(P0+dP1*dt/2,fc0+dfc1*dt/2,fv0+dfv1*dt/2,det,t+dt/2)
    dP3,dfc3,dfv3=sBloch(P0+dP2*dt/2,fc0+dfc2*dt/2,fv0+dfv2*dt/2,det,t+dt/2)
    dP4,dfc4,dfv4=sBloch(P0+dP1*dt/2,fc0+dfc3*dt/2,fv0+dfv3*dt/2,det,t+dt/2)

    P=P0+dt/6*(dP1+dP2+dP3+dP4)
    fc=fc0+dt/6*(dfc1+dfc2+dfc3+dfc4)
    fv=fv0+dt/6*(dfv1+dfv2+dfv3+dfv4)

    return P,fc,fv

N=10000

P=np.zeros(N,dtype=complex)
fc=np.zeros(N,dtype=complex)
fv=np.ones(N,dtype=complex)

t=np.linspace(-3,3,N)

for i in range(1,t.size):
    P[i],fc[i],fv[i]=rk4(P[i-1],fc[i-1],fv[i-1],0,t[i-1],t[1]-t[0])

fig,ax=figure()
ax.plot(t,np.absolute(fc),label=r'electron density')
ax.plot(t,1-np.absolute(fv),label=r'hole density')
ax.plot(t,np.absolute(P),label='Polarization')

ax.legend()
ax.set_xlabel('time (a.u)')
ax.set_ylabel(r'$\rho$ (a.u)')

print(np.absolute(fc+fv))
#print(P)

plt.show()
