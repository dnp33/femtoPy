import numpy as np

def rk4(f,u0,dt,alpha,Imax):
    f1 = f (u0,alpha,Imax)
    f2 = f (u0 + dt * f1 / 2.0,alpha,Imax)
    f3 = f (u0 + dt * f2 / 2.0,alpha,Imax)
    f4 = f (u0 + dt * f3,alpha,Imax)
    d1=( f1 + 2.0 * f2 + 2.0 * f3 + f4 ) / 6.0
    u1 = u0 + dt * d1
    
    return u1,-d1

def sat( I, alpha, Imax):
    return -alpha*Imax*I/(I+Imax)

def absorption(I0, alpha=1, Imax=1, d=np.linspace(0,10,1000)):
    dz=d[1]-d[0]
    I=np.zeros(d.size+1)
    dI=np.zeros(d.size)
    I[0]=I0

    for i in np.arange(1,d.size+1):
        I[i],dI[i-1]=rk4(sat,I[i-1],dz,alpha,Imax)
        # dI[i-1]=I[i-1]-I[i]
    return dI
