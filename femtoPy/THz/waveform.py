from numpy import loadtxt as np_loadtxt; from numpy import amax as np_amax
from numpy import amin as np_amin; from numpy import where as np_where
from numpy import zeros as np_zeros; from numpy import absolute as np_absolute
from numpy import angle as np_angle; from numpy import unwrap as np_unwrap
from numpy import flipud as np_flipud
from numpy.fft import rfft as fft_rfft
from numpy.fft import rfftfreq as fft_rfftfreq
from scipy.signal import hilbert as sgn_hilbert

# class to hold waveform data
class waveform:
    def __init__(self ,folder='.',filename=False,t_col=0,wf_col=1,
                 peak_field=False,tmin=False,tmax=False,tShift=0,
                 fmin=False,fmax=False,t=False,wf=False,dat=False,sens=200,invertTime=False):
        # find data from file
        if filename:
            self.dat=np_loadtxt(folder+'/'+filename)
            self.t=self.dat[:,t_col]+tShift
            
            if invertTime:
                self.t=self.t*-1
                self.wf=self.dat[:,wf_col]/10*sens
            else:
                self.wf=np_flipud(self.dat[:,wf_col])/10*sens
            self.trimTime(tmin,tmax); self.FFT(); self.trimFreq(fmin,fmax)
        # send t & wf data
        if type(wf) != bool:
            self.t=t+tShift; self.wf=wf
            self.trimTime(tmin,tmax); self.FFT(); self.trimFreq(fmin,fmax)
        # send data as single file
        if type(dat) != bool:
            self.dat=dat;
            self.t=-self.dat[:,t_col]+tShift; self.wf=self.dat[:,wf_col]
            self.trimTime(tmin,tmax); self.FFT(); self.trimFreq(fmin,fmax)
        # scale peak_field
        if peak_field:
            self.wf=self.wf/np_amax(self.wf)*peak_field
            self.amp=self.amp/np_amax(self.wf)*peak_field            
        return

    # calculate envelope (does not update automatically after trimming)
    def envelope(self):
        self.env=sgn_hilbert(self.wf)
        self.envPhase=np.angle(self.env)
        self.env=np.absolute(self.env)
        return
    
    # transform
    def FFT(self):
        self.f=fft_rfftfreq(self.t.size,self.t[1]-self.t[0])
        self.amp=fft_rfft(self.wf)
        return
    
    # trim time & waveform array to desired window, calc transform
    def trimTime(self,tmin=False,tmax=False):
        loc_tmin=0; loc_tmax=-1
        if tmin:
            loc_tmin=np_amin(np_where(self.t > tmin))
        if tmax:
            loc_tmax=np_amin(np_where(self.t > tmax))
        if loc_tmax==-1:
            loc_tmax=self.t.size
        self.t[:loc_tmax-loc_tmin]=self.t[loc_tmin:loc_tmax]
        self.t.resize(loc_tmax-loc_tmin,refcheck=False)
        self.wf=self.wf[loc_tmin:loc_tmax]
        self.FFT()
        return
    
    # trim frequencies to useful range
    def trimFreq(self,fmin=False,fmax=False):
        self.loc_fmin=0; self.loc_fmax=-1
        if fmin:
            self.loc_fmin=np_amin(np_where(self.f > fmin))
        if fmax:
            self.loc_fmax=np_amin(np_where(self.f > fmax))
        if self.loc_fmax==-1:
            self.loc_fmax=self.f.size
        self.fTrim=self.f[self.loc_fmin:self.loc_fmax+1]
        self.ampTrim=self.amp[self.loc_fmin:self.loc_fmax+1]
        return
    
    # inverse transform with only useful data range
    def wfTrim(self):
        inv=np_zeros(wlf.f.size); inv[self.loc_fmin:self.loc_fmax]=self.ampTrim
        self.wfTrim=fft.irfft(inv)
        return
    
    # return amplitude & phase
    def Amp(self):
        return np_absolute(self.amp)
    def Phase(self):
        return np_unwrap(np_angle(self.amp))
    def AmpTrim(self):
        return np_absolute(self.ampTrim)
    def PhaseTrim(self):
        return np_unwrap(np_angle(self.ampTrim))
