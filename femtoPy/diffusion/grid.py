import numpy as np

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
    
