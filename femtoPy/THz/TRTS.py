from numpy import sqrt; from numpy import loadtxt
from femtoPy.THz.waveform import waveform
from copy import deepcopy
from importlib import reload

# class to hold data from TRTS
# INPUTS:
# folder/filename = folder/file that has the data
# Col1/Col2 (sens1/sens2) = column (sensitivty) of lock-in 1/2 (500 Hz/250 Hz)
# tmax/tmin/tshift = trim time to tmax/tmin, shift time vector by tshift
# fmin/fmax = min and max frequencies where data is significant
class TRTS:
    def __init__(self,label='data',folder='.',filename=False,Col1=1,Col2=3,
                 sens1=1000,sens2=5,tmin=False,tmax=False,tshift=False,
                 fmin=False,fmax=False,t=False,ref=False,samp=False):
        """test doc_string"""
        if filename:
            self.sens1=sens1; self.sens2=sens2; self.label=label;
            self.Col1=Col1; self.Col2=Col2
        
            # load data, extract reference & sample field, calc FFT
            self.dat=loadtxt(folder+'/'+filename)
            self.t0=-self.dat[:,0]
            L10=self.dat[:,Col1]*sqrt(2)*sens1
            L20=self.dat[:,Col2]*sens2
        
            # create waveform objects (note, extra assignments to save mem.)
            self.ref0=waveform(t=self.t0,wf=L10+L20)
            self.samp0=waveform(t=self.t0,wf=L10-L20)
            self.ref0.t=self.t0; self.samp0.t=self.t0
            self.f0=self.ref0.f; self.samp0.f=self.ref0.f
        if type(t) != bool:
            assert type(ref) != bool, 'must provide reference data'
            assert type(samp) !=bool, 'must provide sample data'
            self.t0=t
            self.ref0=waveform(t=self.t0,wf=ref)
            self.samp0=waveform(t=self.t0,wf=samp)
            self.ref0.t=self.t0; self.samp0.t=self.t0
            self.f0=self.ref0.f; self.samp0.f=self.ref0.f
        
        # trim times & frequencies
        self.ref=deepcopy(self.ref0); self.samp=deepcopy(self.samp0)
        self.trimTime(tmin,tmax)
        self.samp.f=self.ref.f; self.f=self.ref.f; self.trimFreq(fmin,fmax)
        self.calcTrans()
        return
    
    # trim time/freq arrays 
    def trimTime(self,tmin,tmax):
        self.ref.trimTime(tmin,tmax); self.samp.trimTime(tmin,tmax)
        self.t=self.ref.t; self.samp.t=self.t
        return
    def trimFreq(self,fmin,fmax):
        self.ref.trimFreq(fmin,fmax); self.samp.trimFreq(fmin,fmax)
        self.samp.fTrim=self.ref.fTrim; self.fTrim=self.ref.fTrim
        self.calcTrans()
        return
    
    # shift time array for all data
    def tShift(self,tShift):
        self.t[:]=self.t-tShift
        return
    
    # calculate transmission (samp/ref)
    def calcTrans(self):
        self.trans=self.samp.amp/self.ref.amp
        self.transTrim=self.samp.ampTrim/self.ref.ampTrim
        return

