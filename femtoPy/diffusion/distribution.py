'UNITS: micrometer, picosecond, ...'

import numpy as np
import scipy.linalg as linalg
import bandmat as bm


'class to hold the information for each distribution'
class Distribution:
    def __init__(self,grid,weight=1e4,d0=0,args=[],s=8.5e3,u=0.85,A=2.1e-9,B=2.1e-12,C=2.1e-15,T=300,D=0,q=-1):
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
        else: # if function given
            self.density[:,0]=d0(selg.grid.y,args)
        self.i=0
        
        return

    'set up matrices for iteration'
    def prep(self):
        self.LHS=self.grid.I.copy()
        self.RHS=self.grid.I.copy()

    'set up the second order derivative matrix'
    def dif_matrix(self):
        self.LHS=self.LHS.__sub__(self.grid.D2.copy()*self.D[self.i+1])
        self.RHS=self.RHS.__add__(self.grid.D2.copy()*self.D[self.i])
        return

    'Recombination terms (on diagonal matrix elements)'
    def mono(self):
        self.LHS.data[1,:]=self.LHS.data[1,:]+self.grid.dt*self.A/2.
        self.RHS.data[1,:]=self.RHS.data[1,:]-self.grid.dt*self.A/2.
        return
    def bi(self):
        dens=np.asarray(self.dist.density[:,self.i])[:,0]
        dens=dens*dens*self.weight**2
        self.LHS.data[1,:]=self.LHS.data[1,:]+self.grid.dt*dens*self.B/2.
        self.RHS.data[1,:]=self.RHS.data[1,:]-self.grid.dt*dens*self.B/2.
        return
    def tri(self):
        dens=np.asarray(self.dist.density[:,self.i])[:,0]
        dens=dens*dens*dens*self.weight**3
        self.LHS.data[1,:]=self.LHS.data[1,:]+self.grid.dt*dens*self.C/2.
        self.RHS.data[1,:]=self.RHS.data[1,:]-self.grid.dt*dens*self.C/2.
        return
    
    'Source term'
    def source(self,S):
        self.LHS.data[1,:]=self.LHS.data[1,:]-self.grid.dt*S/2.
        self.RHS.data[1,:]=self.RHS.data[1,:]+self.grid.dt*S/2.
        return

    'Electric field term'
    def eField(self,E):
        # Turn E-field into banded matrix
        E_new=self.grid.I.copy()
        E_new.data[1,:]=-np.asarray(E.field[:,self.i])[:,0]
        E_new.data[1,:]=E_new.data[1,:]*self.u*np.sign(self.q)
        # tot. deriv.    M1:           M2:
        # d/dx(E*rho) = E*d/dx rho + (d/dx E)*rho
        M1=self.grid.D1.copy()
        M1=bm.BandMat(1,1,bm.dot_mm(E_new,M1).data[1:4,:])
        M2=self.grid.I.copy()
        M2.data[1,:]=-np.asarray(E.d1[:,self.i])[:,0]
        M2.data[1,:]=M2.data[1,:]*self.grid.dt*self.u*np.sign(self.q)/2.
        self.LHS=self.LHS.__sub__(M1)
        self.LHS=self.LHS.__sub__(M2)
        self.RHS=self.RHS.__add__(M1)
        self.RHS=self.RHS.__add__(M2)
        return

    'Boundary terms'
    def boundary(self):
        C=self.grid.dt*self.D[self.i]/self.grid.dy/self.grid.dy
        self.LHS.data[0,1]=self.LHS.data[0,1]-2*C
        self.LHS.data[1,0]=self.LHS.data[1,0]+2*C*(1+self.grid.dy*self.s/self.D[self.i])
        self.RHS.data[0,1]=self.RHS.data[0,1]+2*C
        self.RHS.data[1,0]=self.RHS.data[1,0]-2*C*(1+self.grid.dy*self.s/self.D[self.i])
        self.LHS.data[1,-1]=2*C
        self.LHS.data[2,-1]=-C
        self.RHS.data[1,-1]=-2*C
        self.RHS.data[2,-1]=C
        return

    def step(self):
        'RHS (matrix)x(vector)'
        self.RHS=bm.dot_mv(self.RHS,np.asarray(self.density[:,self.i])[:,0])
        'solve system'
        self.density[:,self.i+1]=np.asmatrix(linalg.solve_banded((1,1),self.LHS.data,self.RHS,overwrite_ab=True,overwrite_b=True,check_finite=False)).T
        self.i=self.i+1
        
        return

