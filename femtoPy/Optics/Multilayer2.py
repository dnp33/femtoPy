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
    def __init__(self,eps_r=[1],epsTrans=None,epsRefl=None,l=[0],
                 kx=0,ky=0,k0=1):
        """
        multilayer film class

        Notes
        -----
        S12[0,0] != S21[0,0] because the phase shift is different when you are 
        incident on the front vs the back (this makes sense for, e.g, a mirror
        on one side of your domain)

        Parameters
        ----------
        eps_r : complex dilectric function of each layer in a list
        epsTrans : dielectric function in transmission region
        epsRefl : dielectric function in reflection region
        l : thickness of each layer in a list
        kx : x wavevector (for non-normal incidence)
        ky : y wavevector (for non-normal incidence)
        k0 : wavevector (in units of l^-1)

        Returns
        -------
        Sij : 2x2 matrix of 2x2 matrices 
        S11 -> reflection from left
        S22 -> reflection from right
        S21 -> transmission from left to right
        S12 -> transmission from right to left

        device :
           left                                    right
        refl region - layer 1 - layer 2 - ... - trans region

        Sij[0,0]=Px -> Px
           [1,1]=Py -> Py
           [0,1]=Px -> Py
           [1,0]=Py -> Px
        """
        assert type(eps_r)==list,'eps_r must be a list'
        assert type(l)==list,'l must be a list'
        
        self.V_air=-1j*np.matrix([[kx*ky,1+ky**2],[-kx*kx-1,-kx*ky]],
                                 dtype=complex)
        
        # init. and store each layer
        self.layers=np.empty(len(eps_r),dtype=layer)
        for i in range(self.N):
            self.layers[i]=layer(eps_r=eps_r[i],l=l[i],kx=kx,ky=ky,k0=k0,
                                 V_air=self.V_air)
        if type(epsTrans) != type(None):
            self.layers=np.append(self.layers,transRegion(eps_r=epsTrans,
                                              kx=kx,ky=ky,V_air=self.V_air))
        if type(epsRefl) != type(None):
            self.layers=np.append(reflRegion(eps_r=epsRefl,kx=kx,ky=ky,V_air=self.V_air),self.layers)
        # init. global S-matrix
        self.S11=np.matrix([[0,0],[0,0]],dtype=complex)
        self.S12=I.copy()
        self.S22=np.matrix([[0,0],[0,0]],dtype=complex)
        self.S21=I.copy()

        # calculate global scattering matrix
        for i in range(1,self.N):
            self.star_product(i)
        self.star_product2()

        return
    
    def setS(self,S):
        self.S11=S[0]; self.S12=S[1]; self.S21=S[2]; self.S22=S[3]

    # redheffer star product
    def star_product(self,i):
        M1=self.S12*np.linalg.inv((I-self.layers[i].S11*self.S22))
        M2=self.layers[i].S12*np.linalg.inv(I-self.S22*self.layers[i].S11)
        
        self.S11=self.S11+M1*self.layers[i].S11*self.S21
        self.S12=M1*self.layers[i].S12
        self.S21=M2*self.S21
        self.S22=self.layers[i].S11+M2*self.S22*self.layers[i].S12
        
        return
    
    def star_product2(self):
        M1=self.layers[0].S12*np.linalg.inv((I-self.S11*self.layers[0].S22))
        M2=self.S21*np.linalg.inv(I-self.layers[0].S22*self.S11)

        self.S11=self.layers[0].S11+M1*self.S11*self.layers[0].S21
        self.S12=M1*self.S12
        self.S21=M2*self.layers[0].S21
        self.S22=self.S22+M2*self.layers[0].S22*self.S12

        return

    def star_product3(self):
        M1=self.S12*np.linalg.inv((I-self.layers[-1].S11*self.S22))
        M2=self.layers[-1].S21*np.linalg.inv(I-self.S22*self.layers[-1].S11)
        
        self.S11=self.S11+M1*self.layers[-1].S11*self.S21
        self.S12=M1*self.layers[-1].S12
        self.S21=M2*self.S21
        self.S22=self.layers[-1].S22+M2*self.S22*self.layers[-1].S12
        
        return

    @property
    def N(self):
        return len(self.layers)


class transRegion:
    def __init__(self,eps_r=1,kx=0,ky=0,
                 V_air=1j*np.matrix([[0,1],[-1,0]])):
        self.eps_r=eps_r     # dielectric constant
        self.n=np.sqrt(eps_r)    # index of refraction
        self.kz=np.sqrt(eps_r-ky*ky-kx*kx)    # longitudinal wavevector (update here to include magnetic permeability)
        self.kx=kx   # transverse wavevector
        self.ky=ky   # transverse wavevector
        
        # calculate matrices for S-matrix
        self.calcV()
        self.calcS(V_air)

        return
        
    def calcV(self):
        self.V=np.matrix([[self.kx*self.ky,self.eps_r-self.kx**2],
                          [self.ky**2-self.eps_r,-self.kx*self.ky]],
                         dtype=complex)
        self.V=-1j*self.V/self.kz
        self.V_inv=np.linalg.inv(self.V)
        
        return

    def calcS(self,V_air):
        V_air_inv=np.linalg.inv(V_air)
        A=I+V_air_inv*self.V; Ainv=np.linalg.inv(A)
        B=I-V_air_inv*self.V
        #A=I+self.V_inv*self.V; Ainv=np.linalg.inv(A)
        #B=I-self.V_inv*self.V
        
        self.S11=B*Ainv
        self.S12=0.5*(A-B*Ainv*B)
        self.S21=2*Ainv
        self.S22=-Ainv*B

        return

class reflRegion:
    def __init__(self,eps_r=1,kx=0,ky=0,
                 V_air=1j*np.matrix([[0,1],[-1,0]])):
        self.eps_r=eps_r     # dielectric constant
        self.n=np.sqrt(eps_r)    # index of refraction
        self.kz=np.sqrt(eps_r-ky*ky-kx*kx)    # longitudinal wavevector (update here to include magnetic permeability)
        self.kx=kx   # transverse wavevector
        self.ky=ky   # transverse wavevector
        
        # calculate matrices for S-matrix
        self.calcV()
        self.calcS(V_air)

        return
        
    def calcV(self):
        self.V=np.matrix([[self.kx*self.ky,self.eps_r-self.kx**2],
                          [self.ky**2-self.eps_r,-self.kx*self.ky]],
                         dtype=complex)
        self.V=-1j*self.V/self.kz
        self.V_inv=np.linalg.inv(self.V)
        
        return

    def calcS(self,V_air):
        V_air_inv=np.linalg.inv(V_air)
        A=I+V_air_inv*self.V; Ainv=np.linalg.inv(A)
        B=I-V_air_inv*self.V
        #A=I+self.V_inv*self.V; Ainv=np.linalg.inv(A)
        #B=I-self.V_inv*self.V

        self.S11=-Ainv*B
        self.S12=2*Ainv
        self.S21=0.5*(A-B*Ainv*B)
        self.S22=B*Ainv

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
