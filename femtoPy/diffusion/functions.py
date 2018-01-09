import numpy as np
import scipy.linalg as linalg
import time as time  
from scipy.special import erfcx

#####################
# EQUATIONS FOR PDE #
#####################
'diffusion coefficient terms'
def D_d2_dx2(g,mat):
    return vec_times_band(mat.D,g.D2.copy())
def constD_d2_dx2(g,mat):
    return g.D2.copy()*mat.D
def dD_dx_d_dx(g,mat):
    dD_dx=band_times_vec(g.D1.copy(),mat.D)
    return vec_times_band(dD_dx,g.D1.copy())

'Recombination terms'
def A(g,mat):
    return g.I*g.dt*mat.A/2.
def B(g,mat,dist):
    matrix=g.I.copy()
    matrix[1,:]=matrix[1,:]*dist.n()
    return matrix*g.dt*mat.B/2.
def C(g,mat,dist):
    matrix=g.I.copy()
    matrix[1,:]=matrix[1,:]*dist.n()*dist.n()
    return matrix*g.dt*mat.C/2.
def e_h(g,mat,dist1,dist2):
    matrix=g.I.copy()
    matrix[1,:]=matrix[1,:]*dist1.n()*dist2.n()
    return matrix*g.dt*mat.B/2.

'Source term'
def source(g,S):
    return g.I.copy()*S*g.dt/2.

'Electric field terms'
def mu_E_d_dx(g,mat,E):
    return vec_times_band(mat.mu*E,g.D1.copy())
def mu_dE_dx(g,mat,E):
    matrix=g.I.copy()
    matrix[1,:]=matrix[1,:]*mat.mu*band_times_vec(g.D1.copy(),E)
    return matrix
def dmu_dx_E(g,mat,E):
    matrix=g.I.copy()
    matrix[1,:]=matrix[1,:]*band_times_vec(g.D1.copy(),mat.mu)*E
    return matrix

'Boundary terms'
def boundary(G,mat,RHS):
    # t0=time.time()
    C=G.dt/G.dy/G.dy
    C=C*mat.D[0]
    RHS[0,1]=RHS[0,1]+C
    RHS[1,0]=RHS[1,0]-C*(1+G.dy*mat.s/mat.D[0])

    return RHS

def constBoundary(G,mat,RHS):
    # t0=time.time()
    C=G.dt/G.dy/G.dy
    C=C*mat.D

    RHS[0,1]=C
    RHS[1,0]=-C*(1+G.dy*mat.s/mat.D)

    return RHS


###############################
# LINEAR ALGEBRA FOR SOLUTION #
###############################
'RHS matrix times density at current time step'
def calc_RHS(g,matrix,dist):
    dist.i=dist.i+1
    return band_times_vec(g.I.copy()+matrix.copy(),dist.n(step=-1))

'solve linear system'
def solve(g,matrix,vector):
    LHS=g.I.copy()-matrix.copy()
    'solve matrix equation'
    newDens=np.asmatrix(linalg.solve_banded((1,1),LHS,vector,overwrite_ab=False,overwrite_b=False,check_finite=False)).T
    return newDens


###############################
# DIFFUSION COEFFICIENT CALCS #
###############################
'update diffusion coefficient for next step'
def calc_D(dist,mat):
    tau_e=Caughey_Thomas(N=dist.n(),N0=mat.N0,tMax=mat.tMax,tMin=mat.tMin,alpha=mat.alpha)*10**-15
    mat.mu=mobility(tau=tau_e,mstar=mat.mstar*9.11e-31)
    mu_ab=ambipolar_mobility(mu_e=mat.mu)
    mat.D=diffusion_coefficient(mu_ab,dist.T[dist.i])
    mat.D=mat.D*10**3

    
    return 

'Empirical relationship between scattering rate and photoexcited density'
def Caughey_Thomas(N=np.logspace(-4,0,100),N0=1.05,tMax=375,tMin=5,alpha=0.22):
    tau=(tMax-tMin)
    tau=tau/(1+(N/N0)**alpha)+tMin

    return tau

'Carrier mobility as a function of effective mass, charge, and scattering rate'
def mobility(mstar=0.067*9.11e-31,tau=3.2e-13,q=1.6e-19):
    return q*tau/mstar

'calculate diffusion coefficient using mobility and temperature'
def diffusion_coefficient(mu=0.85,T=300):
    return 1.38e-23*T*mu/1.6e-19

'calculate ambipolar mobility from electron and hole mobility'
def ambipolar_mobility(mu_h=0.04,mu_e=0.85):
    return mu_h*mu_e/(mu_h+mu_e)


####################
# RANDOM FUNCTIONS #
####################
'analytic expression for density as a function of time with monomolecular rec.'
def analyticDistribution(x,args):
    #args[0]= time  args[1]=absorption coefficient  
    t=args[0]
    alpha=args[1]
    mat=args[2]
    z=mat.D*t
    
    x1=alpha*np.sqrt(z)-x/(2*np.sqrt(z))
    x2=alpha*np.sqrt(z)+x/(2*np.sqrt(z))
    x3=mat.s*np.sqrt(t/mat.D)+x/2./np.sqrt(z)
    
    
    phi=np.array(0.5*(erfcx(x1)+\
                      (alpha*mat.D+mat.s)/(alpha*mat.D-mat.s)*erfcx(x2))\
                 -mat.s*erfcx(x3)/(alpha*mat.D-mat.s))
    phi[phi==inf]=0

    np.seterr(under='ignore')
    dat=phi*np.exp(-x*x/(4*z)-t/mat.tau)
    np.seterr(under='raise')
    
    dat[np.where(dat==0)]=np.exp(-alpha*x[np.where(dat==0)])
    
    return dat

'banded matrix operations'
def band_times_vec(M,v):
    return np.roll(M[0,:]*v,-1)+M[1,:]*v+np.roll(M[2,:]*v,1)

def vec_times_band(v,M):
    M[0,1:]=M[0,1:]*v[:-1]
    M[1,:]=M[1,:]*v[:]
    M[2,:-1]=M[2,:-1]*v[1:]

    return M



