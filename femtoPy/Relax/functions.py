"""
Module contains a number of commonly used decay functions:

exp : exponential decay
streExp : stretch exponential decay
powLaw : power law decay
"""

from numpy import piecewise as np_piecewise, exp as np_exp

def exp(t,Amp=1,t0=0,tau=1):
    """
    exponential decay for t > t0

    Parameters
    ----------
    t : time array for calc
    Amp : amplitude @ t0, default 1
    t0 : decay start time, default 0
    tau : relaxation time, default 1
    """
    t=t-t0
    return np_piecewise(t,[t < 0, t >=0],
                        [lambda t: 0, lambda t: Amp*np_exp(-t/tau)])

def streExp(t,Amp=1,t0=0,tau=1,beta=1):
    """
    stretch exponential decay for t > t0

    Parameters
    ----------
    t : time array for calc
    Amp : amplitude @ t0, default 1
    t0 : decay start time, default 0
    tau : relaxation time, default 1
    beta : stretch parameter, default 1
    """
    t=t-t0
    return np_piecewise(t,[t < 0, t >=0],
                [lambda t: 0, lambda t: Amp*np_exp(-(t/tau)**beta)])

def powLaw(t,N0=1,b=1,t0=0,c=1):
    """
    power law decay for t > t0

    Parameters
    ----------
    t : time array for calc
    N0 : initial "density", default 1
    b : rate constant, default 1
    t0 : decay start time, default 0
    c : phenomenological power, default 1 for pure bi-molecular)
    """
    t=t-t0
    val=b*t
    N0=N0**(1/c)
    return np_piecewise(val,[val < 0, val >= 0],
                [lambda val: 0, lambda val: 1/(1/N0+val)**c])

def biExp(t,A1=1,tau1=1,A2=1,tau2=10,t0=0):
    """
    bi-exponential decay for t > t0

    Parameters
    ----------
    t : time array for calc
    A1 : amplitude of first exp
    tau1 : decay constant of first exp
    A2 : amplitude of second exp
    tau2 : decay constant of second exp
    t0 : decay start time, default 0
    """
    t=t-t0
    return exp(t,Amp=A1,t0=0,tau=tau1)+exp(t,Amp=A2,t0=0,tau=tau2)
                 
