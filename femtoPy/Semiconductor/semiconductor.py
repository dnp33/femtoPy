import numpy as np

class semiconductor:
    'class which contains coefficients relevant to solve the diffusion equation in semiconductors'
    def __init__(self):
        self.s = 8.5e3
        self.u_ab = 0.0388
        self.D = self.u_ab*1.38e-23*300/1.6e-19
        self.tau = 2.1e-9
        
        return
    def __init__(self,s=8.5e3,u_ab=0.0388,tau=2.1e-9,T=300,D=0):
        self.s=s
        self.u_ab=u_ab
        if type(D) == int:
            if D==0:
                self.D=np.array(u_ab*1.38e-23*T/1.6e-19)
            else:
                self.D=D
        else:
            self.D=D
        self.tau=tau
        
        return]