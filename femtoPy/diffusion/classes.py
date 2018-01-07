'UNITS: micrometer, nanoseconds, ...'
import numpy as np

class distribution:
    def __init__(self,grid,d0=0,T=300):
        self.density=np.matrix(np.empty((grid.depth().size,grid.time().size)))
        
        'density at t=0'
        if type(d0) == int or type(d0) == float:
            self.density[:,0]=np.matrix(np.zeros(grid.y.size)+d0).T
        elif type(d0) == type(self.density):
            self.density[:,0]=d0
        elif type(d0) == type(np.linspace(0,10)):
            self.density[:,0]=np.asmatrix(d0).T

        'set temperature'
        if type(T) == int or type(T) == float:
            self.T=np.zeros(grid.t.size)+T
        else:
            self.T=T

        'initialize step number to 0'
        self.i=0

        return

    def n(self,step=0):
        return np.asarray(self.density[:,self.i+step])[:,0]
    def n_t(self,g,t):
        loc=np.amin(np.where(g.time() >= t))
        return self.density[:,loc]
    def n_z(self,g,z):
        loc=np.amin(np.where(g.depth() >= z))
        return self.density[loc,:]
                    

'class to hold the information for each distribution'
class material:
    def __init__(self,N0=1.05,s=8.5e3,A=1/2.1,B=1./2.1**2,C=1./2.1**3,D=1,q=-1,tMax=375,tMin=5,alpha=0.22,mu=8.5,mstar=0.067):
        'transport parameters'
        self.s=s            # surface recombination velocity
        self.q=q            # charge
        self.D=D            # diffusion coefficient
        self.tMax=tMax      # maximum scattering time
        self.tMin=tMin      # minimum scattering time
        self.alpha=alpha    # constant for caughey_thomas
        self.N0=N0          # constant for caughey_thomas
        self.mu=mu          # mobility
        self.mstar=mstar    # effective mass

        'recombination coefficients'
        self.A=A            # monomolecular
        self.B=B            # bimolecular (or electron/hole)
        self.C=C            # auger
        
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

    def time(self):
        return np.asarray(self.t.T)[:,0]
    
    def depth(self):
        return np.asarray(self.y)[:,0]
