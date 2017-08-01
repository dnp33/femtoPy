from scipy.integrate import cumtrapz
import numpy as np

'class to calculate electric field from 2 charge distributions'
class gauss:
    def __init__(self,grid=gridParams(),E0='undefined'):
        'grid parameters'
        self.grid=grid
        
        'initial field'
        self.field = np.matrix(np.zeros([self.grid.y.size,self.grid.t.size]))
        
        if type(E0) != str:
            if E0.shape==(501,1):
                self.field[:,0]=E0
            elif E0.shape==(1,501):
                self.field[:,0]=E0.T
            else:
                self.field[:,0]=np.asmatrix(E0).T
        
        return

    def step(self,eRho,hRho,i):
        rho=hRho-eRho
        self.field[:,i]=np.asmatrix(cumtrapz(np.append(0,rho),dx=y[1]-y[0],axis=0)).T
        
        return
