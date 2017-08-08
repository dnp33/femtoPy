import numpy as np
import scipy.linalg as linalg
import scipy.sparse as sparse
import bandmat as bm

'class to hold the grid parameters for the FDTD code'
class Grid:
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

'class to hold the information for each distribution'
class Distribution:
    def __init__(self,grid=Grid(),d0=0,args=[],s=8.5e3,u_ab=0.0388,A=2.1e-9,B=2.1e-12,C=2.1e-15,T=300,D=0,q=-1):
        'grid parameters'
        self.grid=grid
        'transport parameters'
        self.s=s         # surface recombination velocity
        self.u_ab=u_ab   # mobility
        self.q=q         # charge
        self.T=T         # Temperature
        if type(D) == int or type(D) == float:
            if D==0:
                self.D=np.array(u_ab*1.38e-23*T/1.6e-19)
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
        self.LHS=self.LHS-self.grid.D2.copy()*self.D[self.i+1]
        self.RHS=self.RHS+self.grid.D2.copy()*self.D[self.i]
        return

    'Recombination terms (on diagonal matrix elements)'
    def mono(self):
        self.LHS[1,:]=self.LHS[1,:]+self.grid.dt*self.A/2.
        self.RHS[1,:]=self.RHS[1,:]-self.grid.dt*self.A/2.
        return
    def bi(self):
        dens=np.asarray(self.dist.density[:,self.i])[:,0]
        dens=dens*dens
        self.LHS[1,:]=self.LHS[1,:]+self.grid.dt*dens*self.B/2.
        self.RHS[1,:]=self.RHS[1,:]-self.grid.dt*dens*self.B/2.
        return
    def tri(self):
        dens=np.asarray(self.dist.density[:,self.i])[:,0]
        dens=dens*dens*dens
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
        M=self.grid.D1.copy()
        M[0,:]=M[0,:]*E*self.q
        M[2,:]=M[2,:]*E*self.q
        self.LHS=self.LHS-M
        self.RHS=self.RHS+M
        return

    'Boundary terms'
    def boundary(self):
        C=self.grid.dt*self.D[self.i]/self.grid.dy/self.grid.dy
        self.LHS[0,1]=self.LHS[0,1]-2*C
        self.LHS[1,0]=self.LHS[1,0]+2*C*(1+self.grid.dy*self.s/self.D[self.i])
        self.RHS[0,1]=self.RHS[0,1]+2*C
        self.RHS[1,0]=self.RHS[1,0]-2*C*(1+self.grid.dy*self.s/self.D[self.i])
        return

    def step(self):
        'RHS (matrix)x(vector)'
        self.RHS=bm.BandMat(1,1,self.RHS)
        self.RHS=bm.dot_mv(self.RHS,np.asarray(self.density[:,self.i])[:,0])
        'solve system'
        self.density[:,self.i+1]=np.asmatrix(linalg.solve_banded((1,1),self.LHS,self.RHS,overwrite_ab=True,overwrite_b=True,check_finite=False)).T
        self.i=self.i+1
        
        return

from scipy.integrate import cumtrapz
'class to calculate electric field from 2 charge distributions'
class Field:
    def __init__(self,grid,E0='undefined'):
        'initial field'
        self.field = np.matrix(np.zeros([grid.y.size,grid.t.size]))
        
        if type(E0) != str:
            if type(E0) == type(self.field):
                self.field[:,0]=E0
            elif type(E0) == int or type(E0) == float:
                self.field[:,0]=self.field[:,0]+E0
            elif E0.size == self.field[:,0].size:
                self.field[:,0]=np.asmatrix(E0)
        return

    def solveGauss(self,grid,rho1,rho2,i,eps_r=1,E_ext='none'):
        rho=rho1.density[:,i]*rho1.q+rho2.density[:,i]*rho2.q/8.85e-12/eps_r
        self.field[:,i]=np.asmatrix(cumtrapz(np.append(0,rho),dx=grid.y[1]-grid.y[0],axis=0)).T
        if E_ext != 'none':
            self.field[:,i]=self.field[:,i]+E_ext
        
        return
