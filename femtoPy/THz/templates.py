from numpy import loadtxt as np_loadtxt,  zeros as np_zeros,\
    linspace as np_linspace, pi as np_pi, sin as np_sin,\
    flipud as np_flipud, sqrt as np_sqrt

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
           nMin=1,nAvg=2,ext='.txt',flip=True):
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

    if flip:
        refWfs=np_flipud(refWfs); sampWfs=np_flipud(sampWfs)

    cls.ref.wfs=refWfs; cls.samp.wfs=sampWfs

    return

def spec1dAvg(cls,file,Col1=1,Col2=3,sens1=1000,sens2=5,nAvg=5,nMin=1):
    """
    import data from L2 lab format

    Parameters
    ----------
    file : path to root file name
    Col1 : lock-in ?? column in data file
    Col2 : lock-in ?? column in data file
    sens1 : r.r/4 ??? lock-in sensitivity 
    sens2 : r.r/2 ??? lock-in sensitivity
    nAvg : number of averages
    nMin : minimum index of average (probably 0 or 1)
    """
    dat=np_loadtxt(file+'_Average.txt')
    cls.ref.t=cls.samp.t=dat[:,0]
    refWfs=np_zeros((cls.t.size,nAvg))
    sampWfs=np_zeros((cls.t.size,nAvg))

    for i in range(nMin,nMin+nAvg):
        dat=np_loadtxt(file+'_run_'+str(i)+'.txt')
        L10=dat[:,Col1]*np_sqrt(2)*sens1
        L20=dat[:,Col2]*sens2
        refWfs[:,i]=L10+L20
        sampWfs[:,i]=L10-L20

    cls.ref.wfs=refWfs; cls.samp.wfs=sampWfs

    return
    
    
def TDS(cls,refFile='',sampFile='',tCol=0,wfCol=1,sensRef=1,sensSamp=1,
         flip=True):
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
    if flip:
        ref[:,wfCol]=np_flipud(ref[:,wfCol])
        samp[:,wfCol]=np_flipud(samp[:,wfCol])
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

def TRTS2dDat(cls,folder):
    """
    This function loads data that has been stored in folder into the 
    spectroscopy2d class

    Notes
    -----
    the naming convention of files in the folder is currently fixed:
    pump time delay -> tPump.dat
    THz-EO/pump delay -> tTHz.dat
    reference waveforms -> ref.dat
    pump waveforms -> pump.dat
    """
    cls.tPump=np_loadtxt(folder+'tPump.dat')
    cls.tTHz=np_loadtxt(folder+'tTHz.dat')
    cls.ref=np_loadtxt(folder+'ref.dat')
    cls.pump=np_loadtxt(folder+'pump.dat')

    return
    
