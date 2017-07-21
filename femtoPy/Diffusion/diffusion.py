import numpy as np
import scipy.sparse as sparse
import scipy.linalg as linalg

'class to hold the FDTD parameters, data, and grid'
class gridParams:
    def __init__(self,dt=0.00003,dy=0.0005,y_min=0,y_max=.1,t_min=0,t_max=.05,d0='undefined',factor=10,args=[]):
        'grid parameters'
        self.dt = dt # grid size for time
        self.dy = dy # grid size for space
        self.y = np.arange(y_min,y_max+dy,dy) # y-array for grid
        self.t = np.arange(t_min,t_max+dt,dt) # t-array for grid
        self.V = np.matrix(np.zeros([self.y.size,self.t.size]))
        
        'Stencil for interior points'
        self.D2_matrix=np.empty((3,self.y.size))
        self.D2_matrix[0,2:]=1
        self.D2_matrix[1,1:-1]=-2
        self.D2_matrix[2,:-2]=1

        'initial density as a function of depth'
        if d0 != 'undefined': # if function given
            self.d0=d0(self.y,args)
        else:  # otherwise start at 0
            self.d0=np.zeros(self.y.size)
        self.V[:,0]=np.matrix(self.d0).T
        
        return
    
    'functions to return max, and min values for the grid'
    def ymax(self):
        return self.np.amax(y)
    def tmax(self):
        return self.np.amax(t)
    def ymin(self):
        return self.np.amin(y)
    def tmin(self):
        return self.np.amin(t)

    'iteratively solve the PDE'
    def solve(self,mat,bulk_rec):
        'iterate the diffusion equation'
        for i in range(0,self.t.size-1):
            # second order derivative matrix (interior points)
            C=mat.D[i]*self.dt/self.dy**2/2.
            matrix=self.D2_matrix*C
            
            # Boundary Conditions
            matrix[0,1]=2*C
            matrix[1,0]=-2*C*(1+self.dy*mat.s/mat.D[i])
            matrix[1,-1]=0
            matrix[2,-1]=0
            
            # Bulk recombination terms
            matrix[1,:]=matrix[1,:]-self.dt*bulk_rec(self.V[:,i],mat)/2.
            
            # LHS
            LHS=-matrix.copy()
            LHS[1,:]=1+LHS[1,:]

            # RHS
            RHS=matrix
            RHS[1,:]=1+RHS[1,:]
            RHS=sparse.diags([RHS[0,1:],RHS[1,:],RHS[2,:-1]],[1,0,-1]).dot(self.V[:,i])
            
            # 
            self.V[:,i+1]=np.asmatrix(linalg.solve_banded((1,1),LHS,RHS,overwrite_ab=True,overwrite_b=True,check_finite=False))

        return
