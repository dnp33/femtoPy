'UNITS: micrometer, nanoseconds, ...'
import numpy as np
import scipy.linalg as linalg
import time as time


'class to hold the information for each distribution'
class Distribution:
    def __init__(self,grid,weight=1,d0=0,args=[],s=8.5e3,u=0.85,A=1/2.1,B=1./2.1**2,C=1./2.1**3,T=300,D=0,q=-1):
        'grid parameters'
        self.grid=grid
        'transport parameters'
        self.weight=weight  # weighted density
        self.s=s            # surface recombination velocity
        self.u=u            # mobility
        self.q=q            # charge
        if type(T) == int or type(T) == float:
            self.T=np.zeros(self.grid.t.size)+T         # Temperature
        else:
            self.T=T
        if type(D) == int or type(D) == float:
            if D==0:
                self.D=u*1.38e-23*T/1.6e-19
            else:
                self.D=D+np.zeros(grid.t.size)
        else:
            self.D=D
        'recombination coefficients'
        self.A=A
        self.B=B
        self.C=C
        'matrix of densities'
        self.density=np.matrix(np.zeros([self.grid.y.size,self.grid.t.size]))

        'initial density as a function of depth'
        if type(d0) == int or type(d0) == float:
            self.density[:,0]=np.zeros(self.grid.y.size)+d0
        elif type(d0) == type(self.density):
            self.density[:,0]=d0
        elif type(d0) == type(np.linspace(0,10)):
            self.density[:,0]=np.asmatrix(d0).T
        else: # if function given
            self.density[:,0]=d0(self.grid.y,args)
        self.i=0

        'TIME FOR EACH EXECUTION'
        self.tprep=0
        self.tdif=0
        self.tboundary=0
        self.tefield=0
        self.tstep=0
        self.test=0
        
        return

    'set up matrices for iteration'
    def prep(self):
        t0=time.time()
        self.LHS=self.grid.I.copy()
        self.RHS=self.grid.I.copy()
        x=time.time()-t0
        self.tprep=self.tprep
        if x==0:
            self.test=self.test+1
        return

    'set up the second order derivative matrix'
    def dif(self):
        t0=time.time()
        self.LHS=self.LHS-self.grid.D2.copy()*self.D[self.i+1]
        self.RHS=self.RHS+self.grid.D2.copy()*self.D[self.i]
        x=time.time()-t0
        self.tdif=self.tdif+x
        if x==0:
            self.test=self.test+1
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
        t0=time.time()

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

        x=time.time()-t0
        self.tefield=self.tefield+x
        if x==0:
            self.test=self.test+1
        
        return

    'Boundary terms'
    def boundary(self,E=0):
        t0=time.time()
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
        x=time.time()-t0
        self.tboundary=self.tboundary+x
        if x==0:
            self.test=self.test+1
        return

    def step(self):
        t0=time.time()
        'RHS (matrix)x(vector)'
        nRHS=self.RHS[1,:]*np.asarray(self.density[:,self.i])[:,0]
        nRHS[:-1]=nRHS[:-1]+self.RHS[0,1:]*np.asarray(self.density[1:,self.i])[:,0]
        nRHS[1:]=nRHS[1:]+self.RHS[2,:-1]*np.asarray(self.density[:-1,self.i])[:,0]

        'solve system'
        self.density[:,self.i+1]=np.asmatrix(linalg.solve_banded((1,1),self.LHS,nRHS,overwrite_ab=True,overwrite_b=True,check_finite=False)).T
        self.i=self.i+1
        x=time.time()-t0
        self.tstep=self.tstep+x
        if x==0:
            self.test=self.test+1
        return

