"""
THz.wvf is a set of classes for waveforms, 1D spectroscopy, and 2D spectroscopy

waveform
--------
single waveform

spectroscopy1d
--------------
waveform + reference

spectroscopy2d
--------------
2d waveform + references

"""
from numpy import amax as np_amax, amin as np_amin, where as np_where,\
    zeros as np_zeros, absolute as np_absolute, angle as np_angle,\
    unwrap as np_unwrap, flipud as np_flipud, sqrt as np_sqrt, \
    mean as np_mean, std as np_std, real as np_real, imag as np_imag,\
    append as np_append, pi as np_pi
from numpy.fft import rfft as fft_rfft, rfftfreq as fft_rfftfreq
from scipy.signal import hilbert as sgn_hilbert
from femtoPy.THz.functions import eps_to_sig as fnc_eps_to_sig,\
    sig_to_eps as fnc_sig_to_eps
from femtoPy.Constants import c 

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
    def __init__(self,loadFunc=None,peak_field=None,sens=None,tmin=None,
                 tmax=None,tShift=0,fmin=None,fmax=None,calc=True,
                 **kwargs):
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
        calc : calculate FFT/wf etc. automatically (default True)
        
        **kwargs : refer to loadFunc
        """
        self.calc=calc; self.avg=False
        
        if type(loadFunc) != type(None):
            loadFunc(self,**kwargs)
            if peak_field: self.scale(peak_field)
            elif sens: self.scale(sens/10)
            else: self.scale(1)
            self.trimTime(tmin=tmin,tmax=tmax)
        
        return

    def shiftTime(self,tShift=None):
        if tShift == 'MAX':
            loc=np_where(np_absolute(self.wf)==np_amax(np_absolute(self.wf)))
            self.t=self.t-self.t[np_amin(loc)]
        elif type(tShift)==type(None): return
        else: self.t=self.t-tShift

        return
        
    def calcWfAvg(self):
        """calculate average of all wfs"""
        self.wf=np_mean(self.wfs,axis=1)
        self.wfErr=np_std(self.wfs,axis=1,ddof=1)/np_sqrt(self.wfs[0,:].size)
        return
    
    def calcFFTavg(self):
        """calculate average of all FFTs"""
        self.fft=np_mean(self.ffts,axis=1)
        self.fftErr=np_std(self.ffts,axis=1,ddof=1)/np_sqrt(self.wfs[0,:].size)
        return
        
    def scale(self,scale):
        """
        scales waveform and wf avgs by scale
        """
        scale=scale/np_amax(np_absolute(self.wf))
        if self.avg:
            self.wfs=self.wfs*scale
            self.calcWfAvg()
        else:
            self.wf=scale*self.wf
        
        if self.calc: self.FFT()

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
    def FFT(self):
        """calculate FFT"""
        self.f=fft_rfftfreq(self.t.size,self.t[1]-self.t[0])
        if self.avg:
            self.ffts=fft_rfft(self.wfs,axis=0)
            self.calcFFTavg()
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
            if self.avg:
                self.wfs=self.wfs[loc,:]
                self.wfErr=self.wfErr[loc]
        if tmax:
            loc=np_where(self.t < tmax)
            self.t=self.t[loc]
            self.wf=self.wf[loc]
            if self.avg:
                self.wfs=self.wfs[loc]
                self.wfErr=self.wfErr[loc]
                                        
        if self.calc: self.FFT()

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
    def t(self):
        return self._t
    @t.setter
    def t(self,t):
        self._t=t
    @property
    def wf(self):
        return self._wf
    @wf.setter
    def wf(self,wf):
        self._wf=wf
        if self.calc: self.FFT()
        return
    @property
    def wfs(self):
        return self._wfs
    @wfs.setter
    def wfs(self,wfs):
        self._wfs=wfs
        self.avg=True
        if self.calc:
            self.calcWfAvg()
            self.FFT()
    @property
    def amp(self):
        """returns fft amplitude"""
        return np_absolute(self.fft)
    @property
    def amps(self):
        """returns amplitudes of each measurement in average"""
        return np_absolute(self.ffts)
    @property
    def ampErr(self):
        """returns error in fft amplitude"""
        return np_absolute(self.fftErr)
    @property
    def phase(self):
        """returns fft phase"""
        return np_unwrap(np_angle(self.fft))
    @property
    def phases(self):
        """returns phase of each measurement in the average"""
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
    def tWfMax(self,peak='absolute'):
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
    def fAmpMax(self):
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
        self.ref=waveform(); self.samp=waveform()
        
        loadFunc(self,**kwargs)
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
        try: self.samp.fftErr
        except: pass
        else:
            self.transErr=self.trans*np_sqrt((self.ref.ampErr/self.ref.amp)**2
                          +(self.samp.ampErr/self.samp.amp)**2)

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
    @property
    def sigma(self):
        """complex conductivity"""
        return self._sigma
    @sigma.setter
    def sigma(self,sigma):
        """
        set complex conductivity & calc epsilon/n
        
        Notes
        -----
        sets epsilon(0) & n(0) to 1 because of the diverging f^-1 term
        """
        self._sigma=sigma
        f,self._eps=fnc_sig_to_eps(self.f,sigma)
        if self._eps.size < self.f.size:
            self._eps=np_append(1,self._eps)
        self.n=np_sqrt(self._eps)
        return
    @property
    def sigmaR(self):
        """real part of conductivity"""
        return np_real(self._sigma)
    @property
    def sigmaI(self):
        """imaginary part of conductivity"""
        return np_imag(self._sigma)
    @property
    def eps(self):
        """complex dielectric function"""
        return self._eps
    @eps.setter
    def eps(self,eps):
        """set dielectric function & calculate conductivity/index"""
        self._eps=eps
        self._sigma=fnc_eps_to_sig(self.f,eps)
        self._n=np_sqrt(eps)
        return
    @property
    def epsR(self):
        """real part of dielectric function"""
        return np_real(self._eps)
    @property
    def epsI(self):
        """imaginary part of dielectric function"""
        return np_imag(self._eps)
    @property
    def n(self):
        """real part of index"""
        return np_real(self._n)
    @n.setter
    def n(self,n):
        self._n=n
        self._eps=n*n
        self._sigma=fnc_eps_to_sig(self.f,self._eps)
        return
    @property
    def k(self):
        """imaginary part of index"""
        return np_imag(np_sqrt(self._n))
    @property
    def nComplex(self):
        """complex index"""
        return self._n
    @property
    def alpha(self):
        """absorption coefficient"""
        return 4*np_pi*self.nI*self.f*1e10/c

    @property
    def sigmaErr(self):
        """error in complex index"""
        return self._sigmaErr
    @sigmaErr.setter
    def sigmaErr(self,sigmaErr):
        """
        set error and calculate index and dielectric function errors
        (calculations incomplete)
        """
        self._sigmaErr=sigmaErr
        return
    @property
    def sigmaErrR(self):
        """error in real index"""
        return np_real(self._sigmaErr)
    @property
    def sigmaErrI(self):
        """error in imaginary index"""
        return np_imag(self._sigmaErr)
    @property
    def epsErr(self):
        """error in complex index"""
        return self._epsErr
    @epsErr.setter
    def epsErr(self,epsErr):
        """
        set dielectric function error and calculate index and conductivity error
        """
        self._epsErr=epsErr
        self._nErr=np.sqrt(epsErr)
        return 
    @property
    def epsErrR(self):
        """error in real part of dielectric function"""
        return np_real(self._epsErr)
    @property
    def epsErrI(self):
        """error in imaginary part of the dielectric function"""
        return np_imag(self._epsErr)
    @property
    def nErr(self):
        """error in complex index"""
        return self._nErr
    @nErr.setter
    def nErr(self,nErr):
        """set index error and calc dielectric function & conductivity error"""
        self._nErr=nErr
        return
    @property
    def nErrR(self):
        """real part of index error"""
        return np_real(self._nErr)
    @property
    def nErrI(self):
        """imaginary part of index error"""
        return np_imag(self._nErr)
