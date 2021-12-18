"""
class for simulating optical transmission through a multilayer film. see Multilayer.multilayer._doc__ for more details on usage.
"""

import numpy as np

I=np.matrix([[1,0],[0,1]],dtype=complex)

class multiLayer:
    def __init__(self,eps_r=[1],epsTrans=None,epsRefl=None,l=[0],
                 kx=0,ky=0,k0=1):
        """
        multilayer film class

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

        typically, we consider light to be incident from left to right so that
	transmission/reflection through a multilayer film is given by:
           left                                    right
        refl region - layer 1 - layer 2 - ... - layer n - trans region
	
	Each Sij is a matrix, the components of which indicate polarizations. 
	Defined with x,y in the plane and z the propagation direction, the following
	convention holds:
        Sij[0,0]=Px -> Px
           [1,1]=Py -> Py
           [0,1]=Px -> Py
           [1,0]=Py -> Px

	** THERE IS NO ANISOTROPY BUILT INTO THE CODE YET, so the off-diagonal
	components should always be zero and the on diagonal components should be 
	the same.

        Usage:
	------
	for transmission from vacuum through a two layer film (eps_1=2+0.1i, l_1=1 
	and eps_2=3+0.2i, l_2=0.5) into a substrate with eps=4 at a wavevector of 1:
	
	import Multilayer as ml
	
	eps_layers=[2+0.1*1j,1+0.2*1j]
	l_layers=[1,0.5]
	eps_refl=1
	eps_trans=4

	k0=1

	film=ml.multilayer(k0=k0,l=l_layers,eps_r=eps_layers,epsRefl=eps_refl
 			,epsTrans=epsTrans)

	t=film.S21[1,1]
	r=film.S11[1,1]

	print('transmissivity: {:.3f}'.format(t))
	print('reflectivity: {:.3f}'.format(r))

        Notes
        -----
        -At one point I was confused because S12[0,0] != S21[0,0] even when the left/right 		dielectric functions were equal. This happens because the phase shift is different
	when you are incident on the front vs the back (this makes sense for, e.g, a mirror
        on one side of your domain)

	-I still have a bit of confusion over the sign convention. I believe that this 
	calculator expects a positive imaginary part of the dielectric constant for 
	absorption and negative imaginary part for gain, but this should be confirmed 
	by the user.
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
        for i in range(1,self.N-1):
            self.setS(self.star_product(self.S,self.layers[i].S))
        self.setS(self.star_product(self.S,self.layers[-1].S))
        self.setS(self.star_product(self.layers[0].S,self.S))

        return

    @property
    def S(self):
        return [[self.S11,self.S12],[self.S21,self.S22]]
    
    def setS(self,S):
        self.S11=S[0][0]; self.S12=S[0][1]; self.S21=S[1][0]; self.S22=S[1][1]
        return
    
    def star_product(self,A,B):
        M1=A[0][1]*np.linalg.inv(I-B[0][0]*A[1][1])
        M2=B[1][0]*np.linalg.inv(I-A[1][1]*B[0][0])

        S11=A[0][0]+M1*B[0][0]*A[1][0]
        S12=M1*B[0][1]
        S21=M2*A[1][0]
        S22=B[1][1]+M2*A[1][1]*B[0][1]
        
        return [[S11,S12],[S21,S22]]

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

    @property
    def S(self):
        return [[self.S11,self.S12],[self.S21,self.S22]]


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
    @property
    def S(self):
        return [[self.S11,self.S12],[self.S21,self.S22]]

        
    
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

    @property
    def S(self):
        return [[self.S11,self.S12],[self.S12,self.S11]]
