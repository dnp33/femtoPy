import numpy as np
import scipy.special as spc
import scipy.integrate as integrate

class PLspec:
    def __init__(self,grid,E=np.linspace(1.35,2.1,100),gamma=0.01,Eg=1.42,
                 dist=0,alpha0=1,theta=1):
        self.T_eV=1.38e-23/1.6e-19
        self.Eg=Eg
        self.E=E
        self.gamma=gamma
        self.alpha0=band_tail_alpha(self.E-self.Eg,gamma=self.gamma,
                                    alpha0=alpha0,theta=theta)
        a_k0=1-np.exp(-self.alpha0)
        self.CONST=E*E*a_k0
        self.dy=grid.dy
        self.PL=np.zeros((E.size,grid.time().size))
        if dist==0:
            self.PL[:,0]=np.zeros(E.size)
        else:
            self.i=0
            self.calcPL(dist)
        self.i=1

        return

    def calcPL(self,dist):
        T=dist.T[dist.i]*self.T_eV
        I0=self.CONST*np.exp(-self.E/T)
        I0=I0/np.sum(I0)
        I=np.zeros(I0.size)
        n=dist.n()
        n=n*n
        for i in range(n.size):
            Inew=-self.alpha0*I+(n[-i]+n[-i-1])*I0
            Inew=Inew*self.dy/2.
            I=(I+Inew)/(1+self.alpha0*self.dy/2.)
        self.PL[:,self.i]=I
        self.i=self.i+1

        return

    def calcPL_NoAbsorption(self,dist):
        T=dist.T[dist.i]*self.T_eV
        I0=self.CONST*np.exp(-self.E/T)
        I0=I0/np.sum(I0)
        I=np.zeros(I0.size)
        n=dist.n()
        n=n*n
        for i in range(n.size):
            Inew=(n[-i]+n[-i-1])*I0
            Inew=Inew*self.dy/2.
            I=(I+Inew)
        self.PL[:,self.i]=I
        self.i=self.i+1

        return

    def calcPL_constT(self,dist):
        T=dist.T*self.T_eV
        I0=self.CONST*np.exp(-self.E/T)
        I0=I0/np.sum(I0)
        I=np.zeros(I0.size)
        n=dist.n()
        n=n*n
        for i in range(n.size):
            Inew=-self.alpha0*I+(n[-i]+n[-i-1])*I0
            Inew=Inew*self.dy/2.
            I=(I+Inew)/(1+self.alpha0*self.dy/2.)
        self.PL[:,self.i]=I
        self.i=self.i+1

        return

    def sum(self):
        self.spec=np.sum(self.PL,axis=1)
        self.int=np.sum(self.spec)

        return

'absorption coefficient with band-tail convolution'
def band_tail_alpha(detuning,gamma=0.005,theta=1,alpha0=1):
    alpha=alpha0*1/(gamma*spc.gamma(1+1./theta))
    
    def f(u,detuning,gamma,theta):
        z=detuning-u
        dat=np.piecewise(z,[ z <= 0, z > 0],[lambda z: 0,lambda z: np.sqrt(z)])
        dat=dat*np.exp(-np.absolute(u/gamma)**theta)
        return dat
    
    dat=np.zeros(detuning.size)
    err=np.zeros(detuning.size)
    np.seterr(under='ignore')
    for i in np.arange(0,detuning.size):
        buf=integrate.quad(f,detuning[i]-gamma*5000,detuning[i],args=(detuning[i],gamma,theta),epsabs=0,epsrel=10**-4)
        dat[i]=buf[0]
        err[i]=buf[1]
    
    alpha=alpha*dat
    
    return alpha

''''' PHOTOLUMINESCENCE '''''
'calculate PL spectrum with a Boltzmann factor for occupancy'
# alpha=absorption coefficient as a function of energy (a.u)
# E=photon energy (a.u)
# T=temperature (must be same units as E)

def boltzmannPL(E,alpha,T=0.025):
    PL=alpha*np.exp(-E/T)
    
    return PL/np.amax(PL)

def band_tail_PL(E,alpha,A=1,mu=0,T=0.024,d=1):
    a_k0=1-np.exp(-alpha*d)
    I=E*E*a_k0
    I=I*np.exp(-(E-mu)/T)
    
    I=I/np.amax(I)
    
    return I*A

'calculate PL spectrum after reabsorption of light leaving a material'
# PL0=normalized PL spectrum
# alpha=absorption coefficient for each component of PL0
# d=grid of positions where PL is calculated
# N=weight of PL emitted at each d
# I=PL intensity
def PL_reabsorption(alpha,PL0,N,dy):
    I=np.zeros(PL0.size)
    for i in range(N.size):
        Inew=-alpha*I+(N[-i]+N[-i-1])*PL0
        Inew=Inew*dy/2.
        I=(I+Inew)/(1+alpha*dy/2.)

    return I
