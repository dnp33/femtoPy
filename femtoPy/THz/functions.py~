from numpy import pi as np_pi; from numpy import conjugate as np_conj
from numpy import where as np_where

epsilon0=8.85e-12
# basic thin film formula
# INPUTS
# n=substrate index (can be imaginary/dispersive)
# d=film thickness (in desired units)
# OUTPUTS
# returns conductivity (in Siemann's/unit(d)
def thinFilm(t,n=2,d=1):
    sigma=(n+1)*(1/t-1)
    return sigma/377/d

# modified thin film formula
# INPUTS
# "" "" 
# sigma=bare conductivity of film
def thinFilm_mod(t,n=2,d=1,sigma=0):
    nEff=n+sigma*d*377
    return thinFilm(t,nEff,d)

# calculate conductivity from Drude model
# INPUTS
# f=frequency
# DC=DC conductivity
# tau=scattering time
def Drude(f,DC=17500,tau=0.18):
    w=2*np_pi*f
    return DC/(1+1j*w*tau)

# calculate dielectric function from Lorentz model
# INPUTS
# epsInf=high frequency dielectric function
# epsSt=low frequency dielectric function
# f0=resonant frequency
# gamma=damping constant
def Lorentz(f,epsInf=1,epsSt=1.1,f0=1.5,gamma=0.1):
    w0=2*np_pi*f0
    w=2*np_pi*f
    return epsInf+(epsSt-epsInf)*w0**2/(w0**2-w**2-1j*gamma*w)

# calculate dielectric function from conductivity
# INPUTS
# f=frequency
# sigma=complex conductivity
def sig_to_eps(f,sigma):
    sigma=np_conjugate(sigma[np_where(f != 0)])
    f=f[np_where(f !=0)]
    w=2*np_pi*f*1e12
    return f,1+1j*sigma/(w*epsilon0)

# calculate conductivity from dielectric function
# INPUTS
# f=frequency
# eps=complex dielectric function
def eps_to_sig(f,eps):
    return -1j*(eps-1)*epsilon0*f*2*np_pi*1e12

# calculate conductivity from Drude model
# INPUTS
# f=frequency
# tau=scattering time
# sigmaDC=DC conductivity
# c=localization parameter
def DrudeSmith(f,tau,sigmaDC,c):
    sigmaDC=(1-c)*sigmaDC
    w=2*np_pi*f
    const=1-1j*w*tau
    return sigmaDC*(1+c/const)/const
