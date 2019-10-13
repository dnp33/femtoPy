from numpy import loadtxt as np_loadtxt,  zeros as np_zeros,\
    linspace as np_linspace, pi as np_pi, sin as np_sin,\
    flipud as np_flipud

def wfDat(cls,t,wf,scale=1):
    """
    this function takes a time & waveform array to store in cls
    
    Parameters
    ----------
    cls : a waveform object
    t : time array
    wf : waveform array
    scale : scales waveform (wf=wf*scale)
    """
    cls.t=t
    cls.wf=wf

    return

def wfAvg(cls,filename='',nMin=1,nAvg=2,
          wfCol=1,tCol=0,err=True):
    """
    this function loads a set of waveforms that will be averaged in cls

    Parameters
    ----------
    cls : waveform object
    filename : root name of files (e.g., ./waveforms/myWf)
    nMin : lowest number to use (probably 0 or 1)
    nAvg : number of averages
    wfCol : coloumn with waveform data
    tCol : column with time data
    """
    dat=np_loadtxt(filename+str(nMin))
    cls.t=dat[:,tCol]
    cls.wfs=np_zeros((cls.t.size,nAvg))

    for i in range(nMin,nMin+nAvg):
        cls.wfs[:,i-nMin]=np_loadtxt(filename+str(i))[:,1]
    
    return

def wf(cls,filename='',wfCol=1,tCol=0):
    """
    Function that loads a single waveform (no averaging)
    
    Parameters
    ----------
    cls : waveform object that will hold the data
    filename : name of waveform file
    wfCol : column that holds waveform data
    tCol : column that holds time data
    """
    dat=np_loadtxt(filename)
    cls.t=dat[:,0]
    cls.wf=dat[:,1]
    
    return

def TDSavg(cls,refFile='',sampFile='',tCol=0,wfCol=1,sensRef=1,sensSamp=1,
           nMin=1,nAvg=2,ext='.txt'):
    """
    this function loads a set of reference & sample waveforms into waveform
    objects that are stored in cls, a spectroscopy1d object

    Parameters
    ----------
    refFile : root filename for reference wf
    sampeFile : " " sample wf
    tCol : time column in wf files
    wfCol : waveform data column in wf files
    sensRef : ref lock-in sensitivity
    sensSamp : sample " "
    nMin : minimum index of average (probably 0 or 1)
    nAvg : total number of averages
    ext : file extension (default .txt)
    """
    dat=np_loadtxt(refFile+str(nMin)+ext)
    cls.ref.t=cls.samp.t=dat[:,tCol]

    refWfs=np_zeros((cls.t.size,nAvg))

    sampWfs=refWfs.copy()
    for i in range(nMin,nMin+nAvg):
        refWfs[:,i-nMin]=np_loadtxt(refFile+str(i)+ext)[:,wfCol]
        sampWfs[:,i-nMin]=np_loadtxt(sampFile+str(i)+ext)[:,wfCol]

    cls.ref.wfs=refWfs; cls.samp.wfs=sampWfs

    return
    
def TDS(cls,refFile='',sampFile='',tCol=0,wfCol=1,sensRef=1,sensSamp=1,
         invertTime=False):
    """
    this function loads a reference & sample waveform into a waveform
    objects that are stored in cls, a spectroscopy1d object

    Parameters
    ----------
    refFile : root filename for reference wf
    sampeFile : " " sample wf
    tCol : time column in wf files (usually 0)
    wfCol : waveform data column in wf files (usually 1)
    sensRef : ref lock-in sensitivity
    sensSamp : sample " "
    invertTime : switch direction of time (necessary for some file types)
    """
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
