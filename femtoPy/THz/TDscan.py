from numpy import empty as np_empty
from numpy import amax as np_amax; from numpy import amin as np_amin
from numpy import absolute as np_abs
from numpy import roll as np_roll
from numpy import where as np_where
from numpy import loadtxt as np_loadtxt
from numpy import meshgrid as np_meshgrid
from numpy.random import random as np_rand
from numpy import arange as np_arange
from scipy.interpolate import interp2d as sp_interp2d
from numpy.fft import rfft as np_rfft
from numpy.fft import irfft as np_irfft
from numpy.fft import rfftfreq as np_rfftfreq

class TDscan:
    def __init__(self,label='data',folder=False,
                 tPump=False,tTHz=False,
                 ref=False,pump=False,tmax=False,tmin=False,
                 fmax=False,fmin=False):
        if folder:
            self.tPump=np_loadtxt(folder+'/tPump.dat')
            self.tTHz=np_loadtxt(folder+'/tTHz.dat')
            self.ref=np_loadtxt(folder+'/ref.dat')
            self.pump=np_loadtxt(folder+'/pump.dat')
            self.avg=(self.ref+self.pump)/2

        if type(tPump)!=bool:
            self.tPump=tPump; self.tTHz=tTHz;

        self.dT=self.pump-self.ref
        
        self.refFFT=np_rfft(self.ref,axis=1)
        self.pumpFFT=np_rfft(self.pump,axis=1)
        self.trans=self.pumpFFT/self.refFFT
        
        self.f=np_rfftfreq(self.tTHz.size,self.tTHz[1]-self.tTHz[0])
        self.fX,self.fY=np_meshgrid(self.f,self.tPump)
        
        return
    
    def calcSigma(self,func,*args,deconv=False):
        if deconv:
            self.sigmaDeconv=np_empty(self.transDeconv.shape,dtype=complex)
            for i in range(self.tPumpDeconv.size):
                self.sigmaDeconv[i,:]=func(self.transDeconv[i,:],*args)
        else:
            self.sigma=np_empty(self.trans.shape,dtype=complex)
            for i in range(self.tPump.size):
                self.sigma[i,:]=func(self.trans[i,:],*args)
        return

    def reGrid(self,noise_dT=3,noise_avg=3):
        self.tPumpSkew,self.dTskew=Skew(self.tTHz,self.tPump,self.dT,noise_dT)
        self.dTSkewFFT=np_rfft(self.dTskew,axis=1)
        
        self.tPumpSkew,self.avgSkew=Skew(self.tTHz,self.tPump,
                                         self.avg,noise_avg)
        self.avgSkewFFT=np_rfft(self.avgSkew,axis=1)
                    
        return

    def deConvolve(self,G_w,noise_dT=3,noise_avg=3,fMax=2.4):
        self.reGrid(noise_dT=noise_dT,noise_avg=noise_avg)
        self.tPumpDeconv=np_arange(np_amin(self.tPump),np_amax(self.tPump),
                                   self.tTHz[1]-self.tTHz[0])
        loc=np_amin(np_where(self.f >= fMax))
        for i in range(self.tPumpSkew.size):
            self.dTSkewFFT[i,:loc]=self.dTSkewFFT[i,:loc]/G_w[:loc]
            self.avgSkewFFT[i,:loc]=self.avgSkewFFT[i,:loc]/G_w[:loc]

        self.dTskew=np_irfft(self.dTSkewFFT,axis=1)
        self.avgSkewFFT=np_irfft(self.avgSkewFFT,axis=1)

        self.dTdeconv=unSkew(self.tTHz,self.tPump,self.tPumpSkew,self.dTskew)
        self.avgDeconv=unSkew(self.tTHz,self.tPump
                              ,self.tPumpSkew,self.avgSkewFFT)
        self.refDeconv=self.avgDeconv-self.dTdeconv
        self.pumpDeconv=self.avgDeconv+self.dTdeconv

        self.refFFTdeconv=np_rfft(self.refDeconv,axis=1)
        self.pumpFFTdeconv=np_rfft(self.pumpDeconv,axis=1)
        self.transDeconv=self.pumpFFTdeconv/self.refFFTdeconv
                                
        return

def unSkew(x,y,ySkew,dat):
    DAT=dat.copy()
    for i in range(x.size):
        DAT[:,i]=np_roll(dat[:,i],i)
        
    return DAT[x.size:,:]

def Skew(x,y,dat,noise=3):
    interp=sp_interp2d(x,y,dat)
    dx=x[1]-x[0]
    ySkew=np_arange(np_amin(y)-dx*x.size,np_amax(y),dx)
    DAT=np_empty((ySkew.size,x.size))

    yMax=np_amax(y); yMin=np_amin(y)
    for i in range(ySkew.size):
        for j in range(x.size):
            if ySkew[i]+j*dx > yMax or ySkew[i]+j*dx < yMin:
                DAT[i,j]=(np_rand(1)-0.5)*noise
            else:
                DAT[i,j]=interp(x[j],ySkew[i]+j*dx)

    return ySkew,DAT
