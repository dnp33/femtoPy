import numpy as np
import numpy.fft as fft

import scipy.optimize as opt

# class to hold data from TRTS
# INPUTS:
# folder/filename = folder/file that has the data
# Col1/Col2 (sens1/sens2) = column (sensitivty) of lock-in 1/2
# tmax/tmin/tshift = trim time to tmax/tmin, shift time vector by tshift
# fmin/fmax = min and max frequencies where data is significant
class TRTS:
    def __init__(self,label='data',folder='',filename='data.dat',Col1=1,Col2=3,sens1=1000,sens2=5,tmin=False,tmax=False,tshift=False,fmin=False,fmax=False):
        self.sens1=sens1; self.sens2=sens2; self.label=label; self.Col1=Col1;
        self.Col2=Col2
        # load data, trim time, extract reference & pumped field
        self.dat=np.loadtxt(folder+'/'+filename)
        self.t0=-self.dat[:,0]
        
        if tmin:
            self.loc_tmin=np.amin(np.where(t0 > tmin))
        else:
            self.loc_tmin=0
        if tmax:
            self.loc_tmax=np.amin(np.where(t0 > tmax))
        else:
            self.loc_tmax=-1
            
        self.L10=dat[:,Col1]*np.sqrt(2)*sens1; self.L20=dat[:,Col2]*sens2
        self.Eref0=L10+L20; self.Epump0=L10-L20

        self.t=self.t[self.loc_tmin:self.loc_tmax]
        self.L1=self.L10[self.loc_tmin:self.loc_tmax]
        self.L2=self.L20[self.loc_tmin:self.loc_tmax]
        self.Eref=L1+L2; self.Epump=L1-L2

        # Calc & trim FFT
        self.f0=fft.rfftfreq(t.size,t[1]-t[0])
        self.Aref0=fft.rfft(self.Eref); self.Apump0=fft.rfft(self.Epump)

        if fmin:
            self.loc_fmin=np.amin(np.where(self.f0 > fmin))
        else:
            self.loc_fmin=0
        if fmax:
            self.loc_fmax=np.amin(np.where(self.f0 > fmax))
        else:
            self.loc_fmax=-1
        self.Aref=self.Aref0[self.loc_fmin:self.loc_fmax]
        self.Apump=self.Apump0[self.loc_fmin:self.loc_fmax]
        self.f=self.f0[self.loc_fmin:self.loc_fmax]

        return 
