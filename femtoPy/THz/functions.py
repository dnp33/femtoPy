"""
THz.functions is set of tools useful for THz spectroscopy

Conductivity Models
-------------------
** 
Note, uses convention that Im{sigma} is positive in the Drude model &
   Im{epsilon} is positive in the Lorentz oscillator model
**
Drude : Drude(f,DC=17500,tau=0.18) returns sigma
Lorentz : Lorentz(f,epsInf=1,epsSt=1.1,f0=1.5,gamma=0.1) returns epsilon
Drude-Smith : DrudeSmith(f,tau,sigmaDC,c)

Conversion Tools
----------------
conductivity to dielectric function : sig_to_eps(f,sigma)
dielectric function to conductivity : eps_to_sig(f,eps)

Spectroscopy Tools
------------------
thinFilm : thinFilm(t,n=2,d=1) n->array or number, d->number
thinFilm_mod :  thinFilm_mod(t,n=2,d=1,sigma=0) n->array or number, d->number, 
                                                sigma->array or number

thinFilmErr : thinFilmErr(sigma,sigmaErr,t,tErr,n,nErr,d,dErr)
"""

from numpy import pi as np_pi, conjugate as np_conj, sqrt as np_sqrt,\
    where as np_where,absolute as np_absolute
from femtoPy.Constants import eps0, imp0 

def thinFilm(trans,n=2,d=1):
    """
    basic thin film formula (needs reference)

    Parameters
    ----------
    n : substrate index (can be imaginary/dispersive)
    d : film thickness (in desired units)
    
    Returns
    -------
    conductivity (in Siemann's/unit(d)
    """
    sigma=(n+1)*(1/trans-1)
    return np_conj(sigma/imp0/d)

def thinFilmErr(sigma,trans,transErr,n,nErr,d,dErr):
    """
    error from thin film formula

    Notes
    -----
    assumes that the conductivity has already been calculated

    Parameters
    ----------
    sigma : conductivity
    trans : transmission
    transErr : error in transmission
    n : substrate index
    nErr : error in substrate index
    d : sample thickness
    dErr : error in sample thickness
    """
    transErr=transErr/(trans**2*(1/trans-1))
    nErr=nErr/(n+1)
    dErr=dErr/d
    sigmaErr=np_absolute(sigma)*np_sqrt(transErr**2+nErr**2+dErr**2)
    
    return sigmaErr

# modified thin film formula
# INPUTS
# "" "" 
# sigma=bare conductivity of film
def thinFilm_mod(t,n=2,d=1,sigma=0):
    """
    modified thin film formula (needs derivation)

    Parameters
    ----------
    n,d : see thinFilm
    sigma : unexcited film conductivity (must be in units(d))
    
    Returns
    -------
    conductivity (in S/unit(d)
    """
    nEff=n+sigma*d*imp0
    return thinFilm(t,nEff,d)

# calculate conductivity from Drude model
# INPUTS
# f=frequency
# DC=DC conductivity
# tau=scattering time
# uses convention that sigma_imag is positive
def Drude(f,DC=17500,tau=0.18):
    """
    Drude conductivity

    Notes
    -----
    uses convention that Im{sigma} is positive
    
    Parameters
    ----------
    f : frequency
    DC : DC conductivity
    tau : scattering time (units(f^-1))
    """
    w=2*np_pi*f
    return DC/(1-1j*w*tau)

# calculate susceptibility from the Lorentz model
# INPUTS
# amp=oscillator strength
# f0=resonant frequency
# gamma=damping constant
def ChiLorentz(f,amp,f0=1.5,gamma=0.1):
    """
    dielectric function of a Lorentz oscillator
    
    Notes
    -----
    uses convention that Im{chi} is positive. The amplitude is currently a.u

    Parameters
    ----------
    amp : oscillator strength
    f0 : resonant frequency
    gamma : damping constant

    Returns
    -------
    epsilon
    """
    w0=2*np_pi*f0
    w=2*np_pi*f
    return amp/(w0**2-w**2-1j*gamma*w)

# calculate dielectric function from Lorentz model
# INPUTS
# epsInf=high frequency dielectric function
# epsSt=low frequency dielectric function
# f0=resonant frequency
# gamma=damping constant
def Lorentz(f,epsInf=1,epsSt=1.1,f0=1.5,gamma=0.1):
    """
    dielectric function of a Lorentz oscillator
    
    Notes
    -----
    uses convention that Im{epsilon} is positive

    Parameters
    ----------
    epsInf : high frequency dielectric function
    epsSt : low frequency dielectric function
    f0 : resonant frequency
    gamma : damping constant

    Returns
    -------
    epsilon
    """
    w0=2*np_pi*f0
    w=2*np_pi*f
    return epsInf+(epsSt-epsInf)*w0**2/(w0**2-w**2-1j*gamma*w)



# INPUTS

def DrudeSmith(f,tau,sigmaDC,c):
    """
    calculate conductivity from Drude-Smith model
    
    Parameters
    ----------
    f : drive frequency
    tau : scattering time
    sigmaDC : DC conductivity
    c : localization parameter

    Returns
    -------
    sigma
    """
    sigmaDC=(1-c)*sigmaDC
    w=2*np_pi*f
    const=1-1j*w*tau
    return sigmaDC*(1+c/const)/const

def sig_to_eps(f,sigma):
    """
    calculate dielectric function from conductivity
    
    Parameters
    ----------
    f : Drive frequency
    sigma : complex conductivity
    
    Returns
    -------
    epsilon
    """
    sigma=sigma[np_where(f != 0)]
    f=f[np_where(f !=0)]
    w=2*np_pi*f*1e12
    return f,1+1j*sigma/(w*eps0)

def eps_to_sig(f,eps):
    """
    calculates the conductivity from the dielectric function

    Parameters
    ----------
    f : drive frequency
    eps : dielectric function
    """
    w=2*np_pi*f*1e12
    return -1j*(eps-1)*eps0*w

# make an FFT figure w/ amplitude & phase
def fftFig(figsize=(10,8),left=.13,right=.87,mid=.7,top=.9,bottom=.13,hspace=0):
    """
    meant to mimic old origin FFT plot style. Mostly unused...
    """
    fig=plt.figure(figsize=figsize)
    'make amplitude figure'
    gs2=gridspec.GridSpec(1,1)
    gs2.update(left=left,right=right,top=mid,bottom=bottom,hspace=0)
    axA=plt.subplot(gs2[0,0])
    axA.set_xlabel('frequency (THz)')
    axA.set_ylabel('amplitude (a.u)')

    'make phase figure'
    gs1=gridspec.GridSpec(1,1)
    gs1.update(left=left,right=right,top=top,bottom=mid)
    axP=plt.subplot(gs1[0,0])
    axP.set_ylabel('phase\n(a.u)')
    # axP.set_xticks([])
    # axP.set_yticks([])

    return fig,axA,axP
