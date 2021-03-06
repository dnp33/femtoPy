# constants
from femtoPy.Constants import ev_nm as c_ev_nm, e as c_e, \
        h as c_h, P as c_P, T as c_T

# misc. ops.
from scipy.optimize import minimize as opt_min
from numpy import loadtxt as np_loadtxt

# math ops.
from numpy import exp as np_exp, average as np_avg,\
        average as np_avg, sqrt as np_sqrt,\
        log as np_log, sum as np_sum
        
# array ops.
from numpy import amax as np_amax, amin as np_amin, \
        where as np_where, zeros as np_zeros

# functions
def Gauss(x,A,x0,w,y0):
    return A*np_exp(-((x-x0)/w)**2)+y0

def loadDPO3032(filename):
    file=open(filename).readlines()

    dat=np_zeros((len(file[16:]),2))
    for i in range(len(file[15:])-1):
        line=file[i+15]
        loc=line.find(',')
        dat[i,0]=line[:loc]
        dat[i,1]=line[loc+1:]
        
    return dat

class Laser:
    def __init__(self,label=''):
        self.label=label
        return
    
    # load spectrum
    def loadSpec(self,filename,bgFile='',bg=False,norm=True):
        dat=np_loadtxt(filename)
        if bg:
            dat[:,1]=dat[:,1]-np_loadtxt(bgFile)[:,1]
            
        locMin=np_amin(np_where(dat[:,0] >= 700))
        locMax=np_amin(np_where(dat[:,0] >= 900))
        dat=dat[locMin:locMax,:]
        
        self.lmda=dat[:,0]; self.E=c_ev_nm/self.lmda
        self.amp=dat[:,1]
        self.lmda0=self.lmda[np_amin(np_where(self.amp==np_amax(self.amp)))]
        self.offset=np_avg(dat[:20],1)
        
        if norm:
            loc=np_amin(np_where(self.amp==np_amax(self.amp)))
            self.A=np_avg(self.amp[loc-8:loc+8])
            self.amp=self.amp/self.A

        return
    
    # load external scope trace
    def loadPDExt(self,HR,LR,pre=0.005,post=0.005):
        HR=loadDPO3032(HR)
        LR=loadDPO3032(LR)
        
        self.tPD=HR[:,0]*1e9; self.LR=LR[:,1]; self.HR=HR[:,1]
        self.tPD=self.tPD-self.tPD[np_amin(np_where(self.LR==np_amax(self.LR)))]
        
        self.PDmax=np_amax(self.LR)
        self.PDpre=pre; self.PDpost=post
        
        self.ratioPre=self.PDmax/self.PDpre
        self.ratioPost=self.PDmax/self.PDpost
        
        return
    # fit spectrum to Gaussian
    def fitSpec(self,p0=[1,800,22,0]):
        # fit function and fit
        def minFunc(p):
            fit=Gauss(self.lmda,p[0],p[1],p[2],p[3])
            return np_sum((fit-self.amp)**2)        
        self.popt=opt_min(minFunc,p0).x
        
        # store values (lmda in nm)
        self.A_fit=self.popt[0]; self.lmda0_fit=self.popt[1]
        self.e=self.popt[2]; self.e2=2*self.popt[2]
        self.offset_fit=self.popt[3]
        self.FWHM=2*self.e*np_sqrt(np_log(2))
        self.fit=Gauss(self.lmda,self.popt[0],self.popt[1],self.popt[2],self.popt[3])
        
        self.E0=c_ev_nm/self.lmda0_fit
        self.f0=self.E0*c_e/c_h/c_T
        
        # calculate energy bandwidth (eV), frequency bandwidth (PHz),
        # and minimum pulse duration (fs)
        self.dE=c_ev_nm/(self.lmda0_fit-self.FWHM/2)-c_ev_nm/(self.lmda0_fit+self.FWHM/2)
        self.df=self.dE*c_e/c_h/c_P
        self.dt_min=0.414/self.df
        
        return
    
    def PrintSpec(self):
        print('--- '+self.label+' ---')
        print('FWHM:\t wavelength\tfrequency\tEnergy')
        print('\t'+str(self.FWHM)[:4],' nm\t',str(self.df*1000)[:5],
              ' THz\t',str(self.dE*1000)[:5]+' meV')
        print('Center:\t wavelength \tfrequency\tEnergy')
        print('\t',str(self.lmda0_fit)[:5],' nm\t',str(self.f0)[:5],
              ' THz\t',str(self.E0)[:3]+' eV')
        print('Transform limit duration - '+str(self.dt_min)[:5]+' fs')
        print('\n')
        return
