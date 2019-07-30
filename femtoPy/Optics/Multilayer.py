"""
class for simulating optical transmission through a multilayer film

INCOMPLETE
"""

import numpy as np

I=np.matrix([[1,0],[0,1]],dtype=complex)
    
# class to hold all the layers in a multilayer structure
# USAGE: INPUTS
# layer (complex) dielectric constant & thickness stored in eps_r & l respectively
# k0 must be in units of l^-1 (e.g., k0=um^-1 & l=um)
# 
# kx & ky in normalized units (i.e., kx=cos(theta)sin(phi), ky=sin(theta)sin(phi))
# set kx=ky=0 for normal incidence (oblique incidence untested)

# USAGE: OUTPUTS
# Sij is a 2x2 matrix
# S11[0,0]=S11[11]=S22[0,0]=S22[1,1]= field transmission coefficient
# S12[0,0]=S12[1,1]= field reflection coefficient
# S21[0,0]=S21[1,1] /= S12[0,0] because the phase shift is different depending whether you're incident from the front or back
class multiLayer:
    def __init__(self,eps_r=[1],l=[0],kx=0,ky=0,k0=1,verbose=False):
        """
        multilayer film class

        Notes
        -----
        S12[0,0] != S21[0,0] because the phase shift is different when you are 
        incident on the front vs the back (I don't remember why...)

        Parameters
        ----------
        eps_r : complex dilectric function of each layer in a list
        l : thickness of each layer in a list
        kx : x wavevector (for non-normal incidence)
        ky : y wavevector (for non-normal incidence)

        Returns
        -------
        Sij -> 2x2 scattering matrix.
        """
        assert type(eps_r)==list,'eps_r must be a list'
        assert type(l)==list,'l must be a list'
        
        self.V_air=-1j*np.matrix([[kx*ky,1+ky**2],[-kx*kx-1,-kx*ky]],dtype=complex)
        
        self.N=len(eps_r)   # number of layers
        
        # init. and store each layer
        self.layers=np.empty(self.N,dtype=layer)
        for i in range(self.N):
            self.layers[i]=layer(eps_r=eps_r[i],l=l[i],kx=kx,ky=ky,k0=k0,V_air=self.V_air) 
        
        # init. global S-matrix (S12 probably equals S21, so this might be more efficient)
        self.S11=np.matrix([[0,0],[0,0]],dtype=complex)
        self.S12=I.copy()
        self.S22=np.matrix([[0,0],[0,0]],dtype=complex)
        self.S21=I.copy()

        # calculate global scattering matrix
        for i in range(self.N):
            self.star_product(i)

        return
    
    # redheffer star product
    def star_product(self,i):
        M1=self.S12*np.linalg.inv((I-self.layers[i].S11*self.S22))
        M2=self.layers[i].S12*np.linalg.inv(I-self.S22*self.layers[i].S11)
        
        self.S11=self.S11+M1*self.layers[i].S11*self.S21
        self.S12=M1*self.layers[i].S12
        self.S21=M2*self.S21
        self.S22=self.layers[i].S11+M2*self.S22*self.layers[i].S12
        
        return
    
# class to hold the information of each layer
# used by multilayer class
class layer:
    def __init__(self,eps_r=1,l=0,kx=0,ky=0,k0=1,V_air=1j*np.matrix([[0,1],[-1,0]])):
        self.eps_r=eps_r     # dielectric constant
        self.n=np.sqrt(eps_r)    # index of refraction
        self.l=l     # layer thickness
        self.kz=np.sqrt(eps_r-ky*ky-kx*kx)    # longitudinal wavevector (update here to include magnetic permeability)
        self.kx=kx   # transverse wavevector
        self.ky=ky   # transverse wavevector
        self.k0=k0   # free space wavevector
        
        # calculate matrices for S-matrix
        self.calcV()
        self.calcX()
        self.calcS(V_air)
        
        return
    
    # update here to include magnetic permeability
    def calcV(self):
        self.V=np.matrix([[self.kx*self.ky,self.eps_r-self.kx**2],[self.ky**2-self.eps_r,-self.kx*self.ky]],dtype=complex)
        self.V=-1j*self.V/self.kz
        self.V_inv=np.linalg.inv(self.V)
        
        return
    
    def calcX(self):
        self.X=I.copy()*np.exp(1j*self.kz*self.k0*self.l)
        
        return
    
    # scattering matrix calculation
    def calcS(self,V_air):
        A=self.V_inv*V_air; B=I-A; A=I+A; Ainv=np.linalg.inv(A)
        M1=self.X*B*Ainv*self.X
        M2=A-M1*B
        M2=np.linalg.inv(M2)
        
        self.S11=M2*(M1*A-B)
        self.S12=M2*self.X*(A-B*Ainv*B)
        
        return
