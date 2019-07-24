from numpy import loadtxt as np_loadtxt,  zeros as np_zeros,\
    linspace as np_linspace, pi as np_pi, sin as np_sin

def wfAvg(cls,filename='',nMin=1,nMax=2,invertTime=False,
          wfCol=1,err=True):
    dat=np_loadtxt(filename+str(1))
    cls.nAvg=nMax-nMin+1
    cls.t=dat[:,0]
    cls.wfs=np_zeros((cls.t.size,cls.nAvg))

    for i in range(nMin,nMax+1):
        cls.wfs[:,i-nMin]=np_loadtxt(filename+str(i))[:,1]
    
    return

def wf(cls,filename):
    dat=np_loadtxt(filename)
    cls.t=dat[:,0]
    cls.wf=dat[:,1]
    
    return

def TDS(cls,refFile,sampFile):
    ref=np_loadtxt(refFile)
    samp=np_loadtxt(sampFile)


def syntheticWaveform(cls,t=np_linspace(-10,10,200),tau=1,t0=0,f0=1,E0=1):
    """
    Function to generate a synthetic pulse for the waveform class
    
    Example
    -------
    wf=THz.waveform.waveform(THz.templates.synthetic)
    
    Parameters
    ----------
    t : time array
    tau : 1/e pulse duration
    t0 : envelope delay
    f0 : center frequency
    E0 : envelope peak
    """
    cls.t=t
    cls.wf=E0*np_exp(-((t-t0)/tau)**2)*np_sin(2*np_pi*f0*t)

    return
    
def syntheticSpec1d(cls):
    """This function is not yet written"""
    return
