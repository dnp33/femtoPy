from numpy import amax as np_amax, amin as np_amin, where as np_where,\
    zeros as np_zeros, absolute as np_absolute, angle as np_angle,\
    unwrap as np_unwrap, flipud as np_flipud, sqrt as np_sqrt,\
    mean as np_mean, std as np_std
from numpy.fft import rfft as fft_rfft, rfftfreq as fft_rfftfreq
from scipy.signal import hilbert as sgn_hilbert

# class to hold waveform data
class waveform:
    def __init__(self,loadFunc,peak_field=None,sens=None,tmin=None,tmax=None,
                 tShift=0,fmin=None,fmax=None,**kwargs):
        try: nAvg
        except: self.nAvg=1
        else: self.nAvg=nAvg
        
        loadFunc(self,**kwargs)
        if peak_field: self.scale(peak_field)
        if sens: self.scale(sens/10)
        if tmin or tmax: self.trimTime(tmin=tmin,tmax=tmax)
        self.t=self.t-tShift

        if self.nAvg > 1:
            self.calcWfAvg()
        
        self.FFT()
        if fmin or fmax: self.trimFreq(fmin=fmin,fmax=fmax)
        
        return

    def calcWfAvg(self,err=True):
        self.wf=np_mean(self.wfs,axis=1)
        if err:
            self.wfErr=np_std(self.wfs,axis=1,ddof=1)/np_sqrt(self.nAvg)
        return
    
    def calcFFTavg(self,err=True):
        self.fft=np_mean(self.ffts,axis=1)
        if err:
            self.fftErr=np_std(self.ffts,axis=1,ddof=1)/np_sqrt(self.nAvg)
        return
        
    def scale(self,scale):
        scale=scale/np_amax(np_absolute(self.wf))
        self.wf=scale*self.wf
        
        try: self.wfs
        except: pass
        else: self.wfs=self.wfs*scale

        return

    # calculate envelope (does not update automatically after trimming)
    def envelope(self):
        self.env=sgn_hilbert(self.wf)
        self.envPhase=np_angle(self.env)
        self.env=np_absolute(self.env)
        return
    
    # transform
    def FFT(self,err=True):
        self.f=fft_rfftfreq(self.t.size,self.t[1]-self.t[0])
        if self.nAvg > 1:
            self.ffts=fft_rfft(self.wfs,axis=0)
            self.calcFFTavg(err=err)
        else:
            self.fft=fft_rfft(self.wf)

        return
    
    # trim time & waveform array to desired window, calc transform
    def trimTime(self,tmin=False,tmax=False):
        if tmin:
            loc=np_where(self.t > tmin)
            self.t=self.t[loc]
            self.wf=self.wf[loc]
            if self.nAvg >1: self.wfs=self.wfs[loc,:]
        if tmax:
            loc=np_where(self.t < tmax)
            self.t=self.t[loc]
            self.wf=self.wf[loc]
            if nAvg > 1: self.wfs=self.wfs[loc]

        return
    
    # trim frequencies to useful range
    def trimFreq(self,fmin=None,fmax=None):
        if fmin:
            loc=np_where(self.f >= fmin)
            self.fTrim=self.f[loc]
            self.fftTrim=self.fft[loc]
        if fmax:
            loc=np_where(self.fTrim <= fmax)
            self.fTrim=self.fTrim[loc]
            self.fftTrim=self.fftTrim[loc]

        return
    
    # return amplitude & phase
    @property
    def amp(self):
        return np_absolute(self.fft)
    @property
    def phase(self):
        return np_unwrap(np_angle(self.fft))
    @property
    def ampTrim(self):
        return np_absolute(self.fftTrim)
    @property
    def phaseTrim(self):
        return np_unwrap(np_angle(self.fftTrim))
