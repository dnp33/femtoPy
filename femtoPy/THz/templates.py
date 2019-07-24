from numpy import loadtxt as np_loadtxt,  zeros as np_zeros

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
