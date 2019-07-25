from numpy import amax as np_amax, amin as np_amin, where as np_where,\
    zeros as np_zeros, absolute as np_absolute, angle as np_angle,\
    unwrap as np_unwrap, flipud as np_flipud, sqrt as np_sqrt,\
    mean as np_mean, std as np_std
from numpy.fft import rfft as fft_rfft, rfftfreq as fft_rfftfreq
from scipy.signal import hilbert as sgn_hilbert

# class to hold waveform data
class waveform:
    """
    This class contains information relevant to time-domain THz waveforms
    & their Fourier transform

    Example
    -------
    wf=THz.waveform.waveform(THz.templates.syntheticWaveform)

    Attributes
    ----------
    **Incomplete section**
    
    """
    def __init__(self,loadFunc=None,peak_field=None,sens=None,tmin=None,tmax=None,
                 tShift=0,fmin=None,fmax=None,**kwargs):
        """
        this __init__ method initializes a waveform object with either a set of
        waveforms from an average or a single waveform. It calculates the 
        Fourier transform, shifts/trims the time vector, and holds an
        additional 'trimmed' frequency array

        Parameters
        ----------
        loadFunc : user defined (or from template module) function of the form 
                   def func(cls, **kwargs):
                        (load waveform data)
                        cls.t=___
                        cls.wf=___  # avg wf. Not required if wfs defined
                        cls.wfs=___ # waveforms from each avergage

        peak_field : scales the waveform to the peak field
        sens : scales waveform to lock-in sensitivity (recessive to peak_field)
        tmin/tmax : min/max time to keep
        fmin/fmax : range of frequencies with good signal/noise
        tShift : time shift
        **kwargs : refer to loadFunc
        """
        try: kwargs['nAvg']
        except: self.nAvg=1
        else: self.nAvg=kwargs['nAvg']
        if type(loadFunc) != type(None):
            loadFunc(self,**kwargs)
            self.calcWfAvg()
            self.shiftTime(tShift)
            if peak_field: self.scale(peak_field,recalc=False)
            elif sens: self.scale(sens/10,recalc=False)
            else: self.scale(1,recalc=False)
            self.trimTime(tmin=tmin,tmax=tmax)
        
        return

    def shiftTime(self,tShift=None):
        if tShift == 'MAX':
            loc=np_where(np_absolute(self.wf)==np_amax(np_absolute(self.wf)))
            self.t=self.t-self.t[np_amin(loc)]
        elif type(tShift)==type(None): return
        else: self.t=self.t-tShift

        return
        
    def calcWfAvg(self,err=True):
        """calculate average of all wfs"""
        self.wf=np_mean(self.wfs,axis=1)
        if err:
            self.wfErr=np_std(self.wfs,axis=1,ddof=1)/np_sqrt(self.nAvg)
        return
    
    def calcFFTavg(self,err=True):
        """calculate average of all FFTs"""
        self.fft=np_mean(self.ffts,axis=1)
        if err:
            self.fftErr=np_std(self.ffts,axis=1,ddof=1)/np_sqrt(self.nAvg)
        return
        
    def scale(self,scale,err=True,recalc=True):
        """
        scales waveform and wf avgs by scale
        """
        scale=scale/np_amax(np_absolute(self.wf))
        
        try: self.wfs
        except: self.wf=scale*self.wf
        else:
            self.wfs=self.wfs*scale
            self.calcWfAvg(err=err)
        
        if recalc:
            self.FFT()

        return

    # calculate envelope (does not update automatically after trimming)
    def calcEnvelope(self):
        """
        Calculates the envelope of the time domain waveform using 
        the Hilbert transform
        """
        self.envlope=sgn_hilbert(self.wf)
        return
    
    # transform
    def FFT(self,err=True):
        """calculate FFT"""
        self.f=fft_rfftfreq(self.t.size,self.t[1]-self.t[0])
        if self.nAvg > 1:
            self.ffts=fft_rfft(self.wfs,axis=0)
            self.calcFFTavg(err=err)
        else:
            self.fft=fft_rfft(self.wf)

        return
    
    # trim time & waveform array to desired window, calc transform
    def trimTime(self,tmin=False,tmax=False):
        """trims time windows, discards old time window"""
        if tmin:
            loc=np_where(self.t > tmin)
            self.t=self.t[loc]
            self.wf=self.wf[loc]
            if self.nAvg >1: self.wfs=self.wfs[loc,:]
        if tmax:
            loc=np_where(self.t < tmax)
            self.t=self.t[loc]
            self.wf=self.wf[loc]
            if self.nAvg > 1: self.wfs=self.wfs[loc]
        self.FFT()

        return

    # trim frequencies to useful range
    def trimFreq(self,fmin=None,fmax=None):
        """creates trimmed frequency & FFT array"""
        if fmin:
            loc=np_where(self.f >= fmin)
            self.fTrim=self.f[loc]
            self.fftTrim=self.fft[loc]
        if fmax:
            loc=np_where(self.fTrim <= fmax)
            self.fTrim=self.fTrim[loc]
            self.fftTrim=self.fftTrim[loc]

        return
    
    @property
    def amp(self):
        """returns fft amplitude"""
        return np_absolute(self.fft)
    @property
    def ampErr(self):
        """returns error in fft amplitude"""
        return np_absolute(self.fftErr)
    @property
    def phase(self):
        """returns fft phase"""
        return np_unwrap(np_angle(self.fft))
    @property
    def ampTrim(self):
        """returns trimmed fft amplitude"""
        return np_absolute(self.fftTrim)
    @property
    def phaseTrim(self):
        """returns trimmed fft phase"""
        return np_unwrap(np_angle(self.fftTrim))
    @property
    def env(self):
        """returns time domain envelope amplitude"""
        return np_absolute(self.envelope)
    @property
    def envPhase(self):
        """returns time domain envelope phase"""
        return np_unwrap(np_angle(self.envelope))
    @property
    def loct_max(self,peak='absolute'):
        """
        returns the location (in time) of the positive/negative/absolute peak
        """
        if peak=='absolute':
            loc=np_where(np_absolute(self.wf)==np_amax(np_absolute(self.wf)))
            return self.t[np_amin(loc)]
        elif peak=='positive':
            return np_amin(np_where(self.wf==np_amax(self.wf)))
        elif peak=='negative':
            return np_amin(np_where(self.wf==np_amin(self.wf)))
        else:
            print('invalid choice for locMax')
            return
    @property
    def locf_max(self):
        """returns the peak frequency"""
        AMP=self.amp
        return f[np_amin(np_where(AMP==np_amax(AMP)))]
    
class spectroscopy1D:
    """
    This class contains two waveform objects: one reference and one sample

    Example
    -------
    spec=THz.wvf.spectroscopy1D(syntheticSpecroscopy)
    """
    def __init__(self,loadFunc,tmin=None,tmax=None,
                 tShift=None,fmin=None,fmax=None,**kwargs):
        """
        this __init__ method initializes two waveform objects (reference &
        sample) to be stored in the 1D-spectroscopy class
        
        Parameters
        ----------
        loadFunc : user defined (or from template module) function of the form
                   def func(ref,samp,**kwargs):
                       (code to load waveform data)
                       ref.t=samp.t=___
                       ref.wf=___ # not required if wfs is defined
                       ref.wfs=___ # 2D array of each average
                       samp.wf=___ # ''
                       samp.wfs=___ # ''
        tmin/tmax,fmin/fmax,tShift : see wvf.waveform
        **kwargs : refer to loadFunc
        """
        try: kwargs['nAvg']
        except: self.nAvg=1
        else: self.nAvg=kwargs['nAvg']
        self.ref=waveform(); self.samp=waveform()
        
        loadFunc(self,**kwargs)
        self.ref.calcWfAvg()
        self.samp.calcWfAvg()
        self.shiftTime(tShift)
        self.trimTime(tmin=tmin,tmax=tmax)
        self.trimFreq(fmin=fmin,fmax=fmax)
        
        return

    def shiftTime(self,tShift):
        """shift time vector of each waveform object (see waveform)"""
        self.ref.shiftTime(tShift)
        self.samp.shiftTime(tShift)
        return

    def trimTime(self,tmin=None,tmax=None):
        """trim time array of each waveform object (see waveform)"""
        self.ref.trimTime(tmin=tmin,tmax=tmax)
        self.samp.trimTime(tmin=tmin,tmax=tmax)

        self.trans=self.samp.fft/self.ref.fft

        return

    def trimFreq(self,fmin=None,fmax=None):
        """trim frequency of each waveform object (see waveform)"""
        self.ref.trimFreq(fmin=fmin,fmax=fmax)
        self.samp.trimFreq(fmin=fmin,fmax=fmax)
        
        if fmin or fmax: self.transTrim=self.samp.fftTrim/self.ref.fftTrim

        return
        
    @property
    def t(self):
        """time vector"""
        return self.ref.t
    @property
    def f(self):
        """frequency vector"""
        return self.ref.f
    @property
    def fTrim(self):
        """trimmed frequency vector"""
        return self.ref.fTrim
