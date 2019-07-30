"tools for Raman spectroscopy. For baseline subtraction, see _A general method for baseline-removal in ultrafast electronpowder diffraction data using the dual-tree complex wavelettransform_, Struc. Dynamics 4, 044004 (2017)"
"code: scikit-ued"
from numpy import loadtxt as np_loadtxt, zeros as np_zeros, \
    append as np_append, flipud as np_flipud, where as np_where,\
    amin as np_amin, amax as np_amax, empty as np_empty, sum as np_sum,\
    ones as np_ones, absolute as np_absolute
from scipy.special import expit as spc_expit
from scipy.sparse import diags as sprs_diags,\
    spdiags as sprs_spdiags
from scipy.sparse.linalg import spsolve as sprs_spsolve
from scipy.optimize import minimize as opt_minimize
from skued import baseline_dt as sku_baseline_dt
from pickle import dump as pk_dump, load as pk_load

def load(filename):
    """
    load a pickled Raman object stored in (filename)
    """
    return pk_load(open(filename+'.fp.raman','rb'))

def baseline_als(y, lam, p, max_iter=10):
    """
    baseline extraction using asymmetric least squares.

    Notes
    -----
    I have no idea how this works but it seems to do a good job. 
    See, e.g., _baseline correction with asymmetric least squares smoothing_, 
    Paul H. C. Eileres & Hans F.M. Boelens
    
    Parameters
    ----------
    y : data with a baseline
    lam : smoothness parameter in some sense... (not really sure)
    p : asymmetry parameter (no idea what this is
    max_iter : number of iterations
    """
    L = len(y)
    D = sprs_diags([1,-2,1],[0,-1,-2], shape=(L,L-2))
    w = np_ones(L)
    for i in range(max_iter):
        W = sprs_spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.transpose())
        z = sprs_spsolve(Z, w*y)
        w = p * (y > z) + (1-p) * (y < z)
    return z

class oscillator:
    """class to hold a Lorentzian function for the Raman class"""
    def __init__(self,wn0,amp,gamm):
        """
        initialize oscillator
        
        Parameters
        ----------
        wn0 : mode wavenumber
        amp : mode amplitude
        gamm : mode linewidth
        """
        self.wn0=wn0
        self.amp=amp
        self.gamm=gamm
        
        return
    def save(self,filename):
        """
        save object to (filename) with python pickle format.
        
        Notes
        -----
        uses extension .fp.ram
        """
        pk_dump(self,open(filename+'.fp.ram','wb'))
        return
    def spec(self,wn):
        """
        calculate the Lorentian lineshape for this mode
        
        Notes
        -----
        the amplitude of the oscillator gives the peak amplitude regardless
        of linewidth. This is not necessarily the correct way to do this...
        
        Parameters
        ----------
        wn : array (or number) over which to calculate the lineshape
        """
        func=1/((wn-self.wn0)**2+self.gamm**2)
        return self.amp*func/np_amax(func)
    
    def label(self,ax,fontsize=15):
        """
        labels the peak corresponding to this mode
        
        Parameters
        ----------
        ax : axis to label
        fontsize : fontisze of label
        """
        return ax.text(self.wn0,self.amp,str(self.wn0),fontsize=fontsize)

class Raman:
    """class to hold tools for Raman spectroscopy"""
    def __init__(self,data,flipud=True):
        """
        initializes a Raman object

        Parameters
        ----------
        data : 2 column array with C0=wn & C1=intensity
        flipud : flip direction of data (e.g., if wn is decreasing with index)
        """
        if flipud: data=np_flipud(data)
        self.wn=data[:,0]; self._I=data[:,1]
        self.baseline=np_zeros(self._I.size)

        self.wn0=np_empty(0); self.amp=np_empty(0);self.gamm=np_empty(0)
        self.modes=np_empty(0,dtype=oscillator)
        
        return
    
    def wvtBase(self,level = 4, max_iter = 30, wavelet = 'qshift1',
                 wnMin=None,wnMax=None,cutoff=0,temp=0.1,baseOffset=0):
        """
        baseline subtraction using dual-tree complex wavelet transform
        
        Parameters
        ----------
        level: decreases "resolution" of baseline, default 4
        max_iter : number of iterations, default 30
        wavelet : choice of wavelet limited by skued_baseline_dt-default qshift1
        wnMin : lower cutoff of wavenumber to include in calc
        wnMax : upper cutoff of wavenumber to include in calc
        cutoff : location of filter cutoff
        baseOffset : baseline in cutoff region
        temp : width of smoothing region
        """
        wn=self.wn; I=self._I
        if type(wnMin) != type(None):
            loc1=np_where(self.wn > wnMin)
            wn=wn[loc1]; I=I[loc1]
            
        if type(wnMax) != type(None):
            loc2=np_where(self.wn < wnMax)
            wn=wn[loc2]; I=I[loc2]

        baseline = sku_baseline_dt(I, level = level, max_iter = max_iter,
                                   wavelet = wavelet)
        baseline=baseline*spc_expit((wn-cutoff)/temp)+\
                  baseOffset*spc_expit(-(wn-cutoff)/temp)
        
        try: loc1
        except: pass
        else: baseline=np_append(np_zeros(np_amin(loc1))+baseOffset,baseline)
        dif=self.wn.size-baseline.size
        if dif > 0: baseline=np_append(baseline,np_zeros(dif))

        self.baseline=baseline

        return

    def alsBase(self,lam=500,p=0.01,max_iter=1000,baseOffset=50,
                cutoff=60,temp=3):
        """
        baseline subtraction using the asymmetric least squares algorithm
    
        Parameters
        ----------
        lam : "order" of filter, default 500 (I don't really understand this)
        p : asymmetry parameter, default 0.01
        max_iter : number of iterations, default 1000
        baseOffset : offset of baseline in filter cuoff region, default 50
        cutoff : location of filter cutoff
        temp : smoothing parameter
        """
        self.baseline=baseline_als(self._I,lam=lam,p=p,max_iter=max_iter)
        self.baseline=self.baseline*spc_expit((self.wn-cutoff)/temp)+\
                  baseOffset*spc_expit(-(self.wn-cutoff)/temp)
        return
    
    def trimWn(self,wnMin=None,wnMax=None):
        """
        trims wavenumber, intensity, & baseline
        
        Parameters
        ----------
        wnMin : minimum desired wavenumber, default None
        wnMax : maximum desired wavenumber, default None
        """
        try: self.baseline
        except: base=False
        else: base=True
        
        if type(wnMin) != type(None):
            loc=np_where(self.wn >=wnMin)
            self.wn=self.wn[loc]; self._I=self._I[loc]
            if base: self.baseline=self.baseline[loc]

        if type(wnMax) != type(None):
            loc=np_where(self.wn <=wnMax)
            self.wn=self.wn[loc]; self._I=self._I[loc]
            if base: self.baseline=self.baseline[loc]

        return
            
    def addMode(self,wn0,amp,gamm):
        """
        adds a mode with wavenumber wn0, amplitude amp, and linewidth gamm
        
        Notes
        -----
        I should add a sorting function
        """
        self.wn0=np_append(self.wn0,wn0)
        self.amp=np_append(self.amp,amp)
        self.gamm=np_append(self.gamm,gamm)
        if type(wn0)==list or type(wn0)==type(self.wn0):
            for i in range(len(wn0)):
                self.modes=np_append(self.modes,oscillator(wn0[i],amp[i],gamm[i]))    
        else: self.modes=np_append(self.modes,oscillator(wn0,amp,gamm))
        return
    def spectrum(self,wn):
        """
        calculates the Raman spectrum resulting from a sum over modes
        
        Parameters
        ----------
        wn : wavenumbers at which to calculate the spectrum (requires array)
        """
        I=np_zeros(wn.size)
        for i in range(self.modes.size):
            I=I+self.modes[i].spec(wn)
            
        return I
    def modModes(self,i=None,wn0=None,amp=None,gamm=None):
        """
        modifies the wavenumber, linewidth, and amplitude of each (or one) mode

        Parameters
        ----------
        i : default None (loops over all modes). Setting i to a number will 
        modify one mode. Can be used as a list to modfiy desired modes
        wavenumber : number (list) with new wavenumber(s)
        amp : " " amp(s)
        gamm : " " gamm(s)
        """
        if type(i)==type(None): i=[i for i in range(self.modes.size)]
        if type(i)==list:
            if type(amp) != type(None):
                for j in i:
                    self.modes[j].amp=amp[j]
                    self.amp[j]=amp[j]
            if type(gamm) != type(None):
                for j in i:
                    self.modes[j].gamm=gamm[j]
                    self.gamm[j]=gamm[j]
            if type(wn0) != type(None):
                for j in i:
                    self.modes[j].wn0=wn0[j]
                    self.wn0[j]=wn0[j]
        elif type(i)==int:
            if type(amp) != type(None):
                self.modes[i].amp=amp
                self.amp[i]=amp
            if type(gamm) != type(None):
                self.modes[i].gamm=gamm[i]
                self.gamm[i]=gamm
            if type(wn0) != type(None):
                self.modes[i].wn0=wn0[i]
                self.wn0[i]=wn0

        return
                
    def label(self,ax,fontsize=15):
        """
        labels an axis with the wavenumber of each mode

        Parameters
        ----------
        ax : axis to label
        fontsize : fontsize of label
        """
        for i in range(self.modes.size):
            Min,Max=ax.get_xlim()
            if self.modes[i].wn0 > Min and self.modes[i].wn0 < Max:
                self.modes[i].label(ax=ax,fontsize=fontsize)
        return

    def fitSpec(self,offset=20,p0=None):
        """
        fit spectrum using the current set of wavenumbers

        Parameters
        ----------
        offset : inital offset guess, default 20
        p0 : new set of amplitudes and linewidths to seed, default None
        """
        def minFunc(p):
            gamm=p[:self.gamm.size]; amp=p[self.amp.size:]
            self.modModes(gamm=gamm,amp=np_absolute(amp))
            fit=self.spec
            
            return np_sum((self.I-offset-fit)**2)
        
        if type(p0)==type(None):
            p0=[i for i in self.gamm]+[j for j in self.amp]
        self.popt=opt_minimize(minFunc,p0)
        self.gamm=self.popt.x[:self.gamm.size]
        self.amp=self.popt.x[self.amp.size:]
        self.offset=offset
        
        return

    def fitWn0(self):
        """
        using the current linewidth & amplitude, optimize the wavenumber
        of each oscillator
        """
        def minFunc(p):
            self.modModes(wn0=p)
            fit=self.spec

            return np_sum((self.I-self.offset-fit)**2)
        p0=self.wn0
        self.wnOpt=opt_minimize(minFunc,p0)
        self.modModes(wn0=self.wnOpt.x)
        
        return
    def printModes(self):
        print('wavenumber\t linewidth\t amplitude')
        for i in range(self.modes.size):
            print('%4.1f'%self.modes[i].wn0,'\t\t','%4.1f'%self.modes[i].gamm,\
                  '\t\t','%4.1f'%self.modes[i].amp)

        return
    @property
    def spec(self):
        """
        returns the calculated Raman spectrum given the current set of modes
        """
        I=np_zeros(self.wn.size)
        for i in range(self.modes.size):
            I=I+self.modes[i].spec(self.wn)

        return I

    @property
    def I(self):
        """
        returns the measured intensity minus the calculated baseline
        """
        return self._I-self.baseline
    
