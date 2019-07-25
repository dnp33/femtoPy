from numpy import loadtxt as np_loadtxt,  zeros as np_zeros,\
    linspace as np_linspace, pi as np_pi, sin as np_sin,\
    flipud as np_flipud

def wfAvg(cls,filename='',nMin=1,nAvg=2,invertTime=False,
          wfCol=1,tCol=0,err=True):
    """
    this function loads a set of waveforms that will be averaged in cls

    Parameters
    ----------
    cls : waveform object
    filename : root name of files (e.g., ./waveforms/myWf)
    nMin : lowest number to use (probably 0 or 1)
    nAvg : number of averages
    invertTime : flip direction of waveform (necessary for some file types)
    wfCol : coloumn with waveform data
    tCol : column with time data
    """
    dat=np_loadtxt(filename+str(nMin))
    cls.nAvg=nAvg
    cls.t=dat[:,tCol]
    cls.wfs=np_zeros((cls.t.size,cls.nAvg))

    for i in range(nMin,nMin+nAvg):
        cls.wfs[:,i-nMin]=np_loadtxt(filename+str(i))[:,1]
    if invertTime: cls.wfs=np_flipud(cls.wfs)
    
    return

def wf(cls,filename='',invertTime=False,wfCol=1,tCol=0):
    """
    Function that loads a single waveform (no averaging)
    
    Parameters
    ----------
    cls : waveform object that will hold the data
    filename : name of waveform file
    invertTime : switch direction of time (necessary for some filetypes)
    wfCol : column that holds waveform data
    tCol : column that holds time data
    """
    dat=np_loadtxt(filename)
    cls.t=dat[:,0]
    cls.wf=dat[:,1]
    if invertTime: cls.wf=np_flipud(cls.wf)
    
    return

def TDSavg(cls,refFile='',sampFile='',tCol=0,wfCol=1,sensRef=1,sensSamp=1,
           nMin=1,nAvg=2,invertTime=False ):
    dat=np_loadtxt(refFile+str(nMin))
    cls.ref.t=cls.samp.t=dat[:,tCol]
    cls.ref.nAvg=nAvg; cls.samp.nAvg=nAvg

    cls.ref.wfs=np_zeros((cls.t.size,cls.nAvg))
    cls.samp.wfs=cls.ref.wfs.copy()
    for i in range(nMin,nMin+nAvg):
        cls.ref.wfs[:,i-nMin]=np_loadtxt(refFile+str(i))[:,wfCol]
        cls.samp.wfs[:,i-nMin]=np_loadtxt(sampFile+str(i))[:,wfCol]

    if invertTime:
        cls.ref.wfs=np_flipud(cls.ref.wfs)
        cls.samp.wfs=np_flipud(cls.samp.wfs)
        
    return
    

def TDS(cls,refFile='',sampFile='',tCol=0,wfCol=1,sensRef=1,sensSamp=1,
         invertTime=False):
    ref=np_loadtxt(refFile)
    samp=np_loadtxt(sampFile)
    cls.ref.t=cls.samp.t=ref[:,tCol]
    cls.ref.wf=ref[:,wfCol]*sensRef
    cls.samp.wf=samp[:,wfCol]*sensSamp

    return


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
