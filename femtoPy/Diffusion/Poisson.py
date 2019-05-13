import numpy as np
from scipy.integrate import cumtrapz
'class to calculate electric field from 2 charge distributions'
class Field:
    def __init__(self,grid,E0='undefined'):
        'initial field'
        self.field = np.matrix(np.zeros([grid.y.size,grid.t.size]))
        self.d1 = np.matrix(np.zeros([grid.y.size,grid.t.size]))
        
        if type(E0) != str:
            if type(E0) == type(self.field):
                self.field[:,0]=E0
            elif type(E0) == int or type(E0) == float:
                self.field[:,0]=self.field[:,0]+E0
            elif E0.size == self.field[:,0].size:
                self.field[:,0]=np.asmatrix(E0)
        return

    def solveGauss(self,rho1,rho2,eps_r=12.9):
        i=rho1.i
        self.d1[:,i]=rho1.density[:,i]*rho1.q+rho2.density[:,i]*rho2.q
        self.d1[:,i]=self.d1[:,i]*rho2.weight/(8.85e-18*eps_r)
        self.field[:,i]=np.asmatrix(cumtrapz(np.append(0,self.d1[:,i]),dx=rho1.grid.dy,axis=0)).T

        return
