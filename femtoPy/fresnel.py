import numpy as np

# Fresnel coefficients
# returns both amplitude and energy transmission/reflection
def rs(n1,n2,thetaI):
    thetaI=thetaI*np.pi/180
    thetaT=snell(n1,n2,thetaI)
    r=(n1*np.cos(thetaI)-n2*np.cos(thetaT))/(n1*np.cos(thetaI)+n2*np.cos(thetaT))
    R=r*r
    return r,R

def rp(n1,n2,thetaI):
    thetaI=thetaI*np.pi/180
    thetaT=snell(n1,n2,thetaI)
    
    r=(n2*np.cos(thetaI)-n1*np.cos(thetaT))/(n2*np.cos(thetaI)+n1*np.cos(thetaT))
    
    return r,r*r

def ts(n1,n2,thetaI):
    thetaI=thetaI*np.pi/180
    thetaT=snell(n1,n2,thetaI)
    t=2*n1*np.cos(thetaI)/(n1*np.cos(thetaI)+n2*np.cos(thetaT))
    T=n2*np.cos(thetaT)*t*t/n1/np.cos(thetaI)
    
    return t,T

def tp(n1,n2,thetaI):
    thetaI=thetaI*np.pi/180
    thetaT=snell(n1,n2,thetaI)
    
    t=2*n1*np.cos(thetaI)/(n2*np.cos(thetaI)+n1*np.cos(thetaT))
    T=n2*np.cos(thetaT)*t*t/n1/np.cos(thetaI)
    return t,T

# calculate transmitted angle using Snell's law
def snell(n1,n2,thetaI):
    return np.arcsin(n1/n2*np.sin(thetaI))

# calculate correction to energy transmission due to the variation in index
def T_cor(n1,n2,thetaI,t):
    thetaT=snell(n1,n2,thetaI)
    return n2*np.cos(thetaT)*t*t/n1/np.cos(thetaI)
