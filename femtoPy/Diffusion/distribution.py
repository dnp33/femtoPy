import scipy.linalg as linalg
import scipy.sparse as sparse

'class to hold the grid parameters for the FDTD code'
class distribution:
    def __init__(self,grid=gridParams(),d0='undefined',args=[],mat=semiconductor()):
        'grid parameters'
        self.grid=grid
        
        'matrix of densities'
        self.density=np.matrix(np.zeros([grid.y.size,grid.t.size]))
        
        'carrier parameters'
        self.mat=mat
        
        'Derivative stencils for interior points'
        self.D1_matrix=np.zeros((3,self.grid.y.size))
        self.D1_matrix[0,2:]=1
        self.D1_matrix[2,:-2]=-1
        self.D1_matrix=self.D1_matrix/4./self.grid.dy
        
        self.D2_matrix=np.zeros((3,self.grid.y.size))
        self.D2_matrix[0,2:]=1
        self.D2_matrix[1,1:-1]=-2
        self.D2_matrix[2,:-2]=1
        self.D2_matrix=self.D2_matrix*self.grid.dt/self.grid.dy**2/2.
        
        'initial density as a function of depth'
        if d0 != 'undefined': # if function given
            self.d0=d0(self.grid.y,args)
        else:  # otherwise start at 0
            self.d0=np.zeros(self.grid.y.size)
        self.density[:,0]=np.matrix(self.d0).T
        
        return

    'iteratively solve the PDE'
    def step(self,bulk_rec,i,E_field):
        'step forward diffusion equation in time'
        # second derivative matrix (interior points)
        C=self.mat.D[i]*self.grid.dt/self.grid.dy**2/2.
        matrix=self.D2_matrix*self.mat.D[i]
        
        # first derivative matrix (interior points)
        buf=self.D1_matrix
        buf[0,:]=buf[0,:]*E_field*self.mat.q
        buf[2,:]=buf[2,:]*E_field*self.mat.q
        matrix=matrix-buf

        # Boundary Conditions
        matrix[0,1]=2*C
        matrix[1,0]=-2*C*(1+self.grid.dy*self.mat.s/self.mat.D[i])
        matrix[1,-1]=0
        matrix[2,-1]=0

        # Bulk recombination and current terms
        matrix[1,:]=matrix[1,:]-self.grid.dt*bulk_rec(self.density[:,i],self.mat)/2.

        # LHS
        LHS=matrix.copy()
        LHS[1,:]=1-LHS[1,:]

        # RHS
        RHS=matrix
        RHS[1,:]=1+RHS[1,:]
        RHS=sparse.diags([RHS[0,1:],RHS[1,:],RHS[2,:-1]],[1,0,-1]).dot(self.density[:,i])

        # 
        self.density[:,i+1]=np.asmatrix(linalg.solve_banded((1,1),LHS,RHS,overwrite_ab=True,overwrite_b=True,check_finite=False))

        return
