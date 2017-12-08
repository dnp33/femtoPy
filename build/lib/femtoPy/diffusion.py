'UNITS: micrometer, nanoseconds, ...'
import numpy as np
import scipy.linalg as linalg
import time as time
from scipy.special import erfcx


'class to hold the information for each distribution'
class dist:
    def __init__(self,grid,weight=1,d0=1,s=8.5e3,A=1/2.1,B=1./2.1**2,C=1./2.1**3,D=1,q=-1):
        'grid parameters'
        self.grid=grid
        'transport parameters'
        self.weight=weight  # weighted density
        self.s=s            # surface recombination velocity
        self.q=q            # charge
        if type(D) == int or type(D) == float:
            self.D=D+np.zeros(grid.t.size)
        else:
            self.D=D
        'recombination coefficients'
        self.A=A
        self.B=B
        self.C=C
        'matrix for density as a function of grid position'
        self.density=np.matrix(np.zeros([self.grid.y.size,self.grid.t.size]))

        'initial condition for density'
        if type(d0) == int or type(d0) == float:
            self.density[:,0]=np.zeros(self.grid.y.size)+d0
        elif type(d0) == type(self.density):
            self.density[:,0]=d0
        elif type(d0) == type(np.linspace(0,10)):
            self.density[:,0]=np.asmatrix(d0).T
        self.i=0

        # 'TIME FOR EACH EXECUTION'
        # self.tdif=0
        # self.tboundary=0
        # self.tefield=0
        # self.tstep=0
        
        return

    'set up matrices for iteration'
    def prep(self):
        # t0=time.time()
        self.LHS=self.grid.I.copy()
        self.RHS=self.grid.I.copy()
        # for timing, commented out..
        # x=time.time()-t0
        # self.tprep=self.tprep
        # if x==0:
        #     self.test=self.test+1
        return

    'set up the second order derivative matrix'
    def dif(self):
        # t0=time.time()
        # M1=self.grid.D2.copy()
        self.LHS=self.LHS-self.grid.D2.copy()*self.D[self.i+1]
        self.RHS=self.RHS+self.grid.D2.copy()*self.D[self.i]
        # x=time.time()-t0
        # self.tdif=self.tdif+x
        # if x==0:
        #     self.test=self.test+1
        return

    'Recombination terms (on diagonal matrix elements)'
    def mono(self):
        self.LHS[1,:]=self.LHS[1,:]+self.grid.dt*self.A/2.
        self.RHS[1,:]=self.RHS[1,:]-self.grid.dt*self.A/2.
        return
    def bi(self):
        dens=np.asarray(self.density[:,self.i])[:,0]
        self.LHS[1,:]=self.LHS[1,:]+self.grid.dt*dens*self.B/2.
        self.RHS[1,:]=self.RHS[1,:]-self.grid.dt*dens*self.B/2.
        return
    def tri(self):
        dens=np.asarray(self.density[:,self.i])[:,0]
        dens=dens*dens*self.weight**3
        self.LHS[1,:]=self.LHS[1,:]+self.grid.dt*dens*self.C/2.
        self.RHS[1,:]=self.RHS[1,:]-self.grid.dt*dens*self.C/2.
        return
    
    'Source term'
    def source(self,S):
        self.LHS[1,:]=self.LHS[1,:]-self.grid.dt*S/2.
        self.RHS[1,:]=self.RHS[1,:]+self.grid.dt*S/2.
        return

    'Electric field term'
    def eField(self,E):
        # t0=time.time()

        E_mu_q=-np.asarray(E.field[:,self.i])[:,0]*self.u*np.sign(self.q)

        # tot. deriv.    M1:           M2:
        # d/dx(E*rho) = E*d/dx rho + (d/dx E)*rho
        # THESE EQUATIONS SHOULD BE CHECKED FOR ACCURACY.....
        M1=self.grid.D1.copy()
        M1[0,:]=M1[0,1:]*E_mu_q[:-1]
        M1[2,:]=M1[2,:-1]*E_mu_q[1:]

        M2=self.grid.D1.copy()
        M2[0,:]=M2[0,1:]*E_mu_q[1:]
        M2[2,:]=M2[2,:-1]*E_mu_q[:-1]

        self.LHS=self.LHS-M1-M2
        self.RHS=self.RHS+M1+M2

        # x=time.time()-t0
        # self.tefield=self.tefield+x
        # if x==0:
        #     self.test=self.test+1
        
        return

    'Boundary terms'
    # CURRENTLY TREATS SECOND BOUNDARY AS ALWAYS 0 
    def boundary(self,E=0):
        # t0=time.time()
        C=self.grid.dt/self.grid.dy/self.grid.dy
        CL=C*self.D[self.i+1]
        CR=C*self.D[self.i]
        self.LHS[0,1]=self.LHS[0,1]-CL
        self.LHS[1,0]=self.LHS[1,0]+CL*(1+self.grid.dy*self.s/self.D[self.i+1])
        self.RHS[0,1]=self.RHS[0,1]+CR
        self.RHS[1,0]=self.RHS[1,0]-CR*(1+self.grid.dy*self.s/self.D[self.i])
        
        # 'sec. boudnary treated the same as the other side w/o surf. rec. term'
        # self.LHS.data[1,-1]=2*C
        # self.LHS.data[2,-1]=-C
        # self.RHS.data[1,-1]=-2*C
        # self.RHS.data[2,-1]=C
        # 'curvature at sec. boundary equal to curvature just inside boundary'
        # buf=bm.BandMat(2,1,np.zeros([4,self.grid.y.size]))
        # buf.data[3,-1]=CL/2.
        # buf.data[2,-1]=-CL
        # buf.data[1,-1]=CL/2.
        # self.LHS=self.LHS.__add__(buf)
        # self.RHS=self.RHS.__sub__(buf)
        'NEITHER OF THESE WORKED........'

        # timing, commented out..
        # x=time.time()-t0
        # self.tboundary=self.tboundary+x
        # if x==0:
        #     self.test=self.test+1
        return

    def step(self):
        # t0=time.time()
        'RHS (matrix)x(vector)'
        nRHS=self.RHS[1,:]*np.asarray(self.density[:,self.i])[:,0]
        nRHS[:-1]=nRHS[:-1]+self.RHS[0,1:]*np.asarray(self.density[1:,self.i])[:,0]
        nRHS[1:]=nRHS[1:]+self.RHS[2,:-1]*np.asarray(self.density[:-1,self.i])[:,0]

        'solve matrix equation'
        self.density[:,self.i+1]=np.asmatrix(linalg.solve_banded((1,1),self.LHS,nRHS,overwrite_ab=True,overwrite_b=True,check_finite=False)).T
        self.i=self.i+1

        # timing, commented out
        # x=time.time()-t0
        # self.tstep=self.tstep+x
        # if x==0:
        #     self.test=self.test+1
        return

'class to hold the grid parameters for FDTD simulation'
class grid:
    def __init__(self,dt=0.05,dy=0.01,y_min=0,y_max=5,t_min=0,t_max=100):
        'grid parameters'
        self.dt = dt # grid size for time
        self.dy = dy # grid size for space
        self.y = np.asmatrix(np.arange(y_min,y_max+dy,dy)).T # y-array for grid
        self.t = np.asmatrix(np.arange(t_min,t_max+dt,dt)) # t-array for grid

        'Derivative stencils for interior points'
        self.D1=np.zeros((3,self.y.size))
        self.D1[0,2:]=1
        self.D1[2,:-2]=-1
        self.D1=self.D1*self.dt/4./self.dy
        
        self.D2=np.zeros((3,self.y.size))
        self.D2[0,2:]=1
        self.D2[1,1:-1]=-2
        self.D2[2,:-2]=1
        self.D2=self.D2*self.dt/self.dy**2/2.

        'Identity matrix'
        self.I=np.zeros((3,self.y.size))
        self.I[1,:]=self.I[1,:]+1
        
        return

'Empirical relationship between scattering rate and photoexcited density'
def Caughey_Thomas(N=np.logspace(14,18,100),N0=1.05e17,tMax=375,tMin=5,alpha=0.22):
    tau=(tMax-tMin)
    tau=tau/(1+(N/N0)**alpha)+tMin

    return tau

'Carrier mobility as a function of effective mass, charge, and scattering rate'
def mobility(mStar=0.067*9.11e-31,tau=3.2e-13,q=1.6e-19):
    return q*tau/mStar

def diffusion_coefficient(mu=0.85,T=300):
    return 1.38e-23*T*mu/1.6e-19

def ambipolar_mobility(mu_h=0.04,mu_e=0.85):
    return mu_h*mu_e/(mu_h+mu_e)

'analytic expression for density as a function of time with monomolecular rec.'
def analyticDistribution(x,args):
    #args[0]= time  args[1]=absorption coefficient  
    t=args[0]
    alpha=args[1]
    mat=args[2]
    z=mat.D*t
    
    x1=alpha*np.sqrt(z)-x/(2*np.sqrt(z))
    x2=alpha*np.sqrt(z)+x/(2*np.sqrt(z))
    x3=mat.s*np.sqrt(t/mat.D)+x/2./np.sqrt(z)
    
    
    phi=np.array(0.5*(erfcx(x1)+\
                      (alpha*mat.D+mat.s)/(alpha*mat.D-mat.s)*erfcx(x2))\
                 -mat.s*erfcx(x3)/(alpha*mat.D-mat.s))
    phi[phi==inf]=0

    np.seterr(under='ignore')
    dat=phi*np.exp(-x*x/(4*z)-t/mat.tau)
    np.seterr(under='raise')
    
    dat[np.where(dat==0)]=np.exp(-alpha*x[np.where(dat==0)])
    
    return dat
