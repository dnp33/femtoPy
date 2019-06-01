import numpy as np
import numpy.fft as fft

import scipy.optimize as opt

class wf:
    'load and transform data'
    def __init__(self,refFile='ref.txt',sampFile='samp.txt',tFile='t.txt',delta=1,l=0.4,c=0.299796):
        'load data files'
        self.refT0=np.loadtxt(refFile)
        self.sampT0=np.loadtxt(sampFile)
        self.sampT=np.copy(self.sampT0)
        self.t=np.loadtxt(tFile)
        
        'transforms'
        self.f0=fft.rfftfreq(self.t.size,self.t[1]-self.t[0])
        self.w0=2*np.pi*self.f0
        self.refF0=fft.rfft(self.refT0)
        self.sampF0=fft.rfft(self.sampT0)
        
        self.f=self.f0
        self.w=self.w0
        self.refF=self.refF0
        self.sampF=self.sampF0
        
        'experimental transfer function'
        self.transE=self.sampF0/self.refF0
        
        'constants'
        self.nAir=np.ones(self.f.size)
        self.l=l
        self.delta=delta
        self.c=c
        
        return
    
    def setValues(self,l='undefined',delta='undefined',nAir='undefined'):
        if nAir != 'undefined':
            self.nAir=np.zeros(self.f.size)+nAir
        if l != 'undefined':
            self.l=l
        if delta != 'undefined':
            self.delta=delta
            
        return
    
    'trim time data and retransform'
    def trimTime(self,tMin=-10,tMax=25,delta='undefined'):
        if tMin > np.amin(self.t):
            loc1=np.amin(np.where(self.t > tMin))
        else:
            loc1=0
        if tMax < np.amax(self.t):
            loc2=np.amin(np.where(self.t > tMax))
        else:
            loc2=self.t.size
            
        self.refT0=self.refT0[loc1:loc2]
        self.sampT0=self.sampT0[loc1:loc2]
        self.sampT=np.copy(self.sampT0)
        self.t=self.t[loc1:loc2]
        
        self.f0=fft.rfftfreq(self.t.size,self.t[1]-self.t[0])
        self.w0=2*np.pi*self.f0
        self.refF0=fft.rfft(self.refT0)
        self.sampF0=fft.rfft(self.sampT0)
        
        self.f=self.f0
        self.w=self.w0
        self.refF=self.refF0
        self.sampF=self.sampF0
        
        self.transExp=self.sampF0/self.refF0
        
        self.nAir=np.ones(self.f.size)
        
        if delta != 'undefined':
            self.delta=delta
        
        return

    'trim frequency and transform arrays'
    def trimFreq(self,fMin=0,fMax='undefined'):
        if fMax=='undefined':
            self.fMax=np.amax(self.f0)
        else:
            self.fMax=fMax
            
        self.fMin=fMin

        self.locMin=np.amin(np.where(self.f0 > self.fMin))
        self.locMax=np.amin(np.where(self.f >= self.fMax))
        
        self.f=self.f0[self.locMin:self.locMax]
        self.w=2*np.pi*self.f
        self.refF=self.refF0[self.locMin:self.locMax]
        self.sampF=self.sampF0[self.locMin:self.locMax]

        self.sampT=np.copy(self.sampF0)
        self.sampT[:self.locMin]=0
        self.sampT[self.locMax:]=0
        
        self.sampT=fft.irfft(self.sampT)
        
        self.nAir=np.ones(self.f.size)
        self.sampCalc=np.copy(self.refF0)        
        
        return
    
    def fPeak(self):
        return self.f[np.amin(np.where(np.absolute(self.refF)==np.amax(np.absolute(self.refF))))]
        
    def calcSamp(self,n,k):
        self.transC=trans(self.nAir,n-1j*k,self.w,self.l,self.delta,self.c)
        self.sampCalc[self.locMin:self.locMax]=self.transC*self.refF
        self.sampTC=fft.irfft(self.sampCalc)
        
        return
    
    def dif(self):    
        overlap=np.sum(np.absolute(self.sampT-self.sampTC))
        integral=np.absolute(np.sum(np.absolute(self.sampT))-np.sum(np.absolute(self.sampTC)))
        
        return overlap+integral

def transAnal(n1,n2,w,l,c):
    return T(n1,n2)*T(n2,n1)*P(n2,l,w,c)*P(n1,-l,w,c)/(1-P(n2,l,w,c)**2*R(n2,n1)**2)

def unwrap(f,phase,loc):
    newPhase=np.zeros(phase.size)
    newPhase[loc]=phase[loc]
    n=0
    if loc < phase.size-1:
        for i in range(loc+1,phase.size):
            if np.absolute(phase[i-1]-phase[i]) > np.pi:
                newPhase[i:]=newPhase[i:]+2*np.pi*np.sign(newPhase[i-1]-newPhase[i-2])
#                 ax.axvline(x[i],color='k',linewidth=0.5,linestyle='--')
            newPhase[i]=phase[i]+newPhase[i]
    if loc > 0:
        for i in range(loc-1,-1,-1):
            if np.absolute(phase[i+1]-phase[i]) > np.pi:
                newPhase[:i+1]=newPhase[:i+1]-2*np.pi*np.sign(newPhase[i+2]-newPhase[i+1])
#                 ax.axvline(x[i],color='k',linewidth=0.5,linestyle='--')
            newPhase[i]=phase[i]+newPhase[i]
    
    n=10
    intercept=np.polyfit(f[:n],newPhase[:n],1)[1]
    newPhase=newPhase-intercept
    
    return newPhase

# amplitude reflection, transmission, and propagation coefficients for normal incidence
def R(n1,n2):
    return (n2-n1)/(n2+n1)

def T(n1,n2):
    return 2*n1/(n1+n2)

def P(n,l,w,c):
    return np.exp(-1j*w*n*l/c)

# n1: frequency dependent index of first medium (probably air)
# n2: frequency dependent index of substrate
# w: angular frequency
# l: material thickness
# numIt: number of iterations to keep in loop
def trans(n1,n2,w,l,numIt,c):
    amp=np.ones(w.size,dtype=complex)
    R2=R(n2,n1)**2
    P2=P(n2,l,w,c)**2
    
    for i in range(1,numIt):
        amp=amp+(P2*R2)**i
    return amp*T(n1,n2)*T(n2,n1)*P(n2,l,w,c)*P(n1,-l,w,c)
