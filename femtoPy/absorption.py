import numpy as np

''''' ABSORPTION OF LIGHT '''''
'saturation model for absorption of light'
# alpha = absorption coefficient at 0 intensity (a.u)
# I = intensity of light (a.u)
# Imax = maximum intensity of light absorbed at a given depth (same units as I)
def satModel( I, alpha, Imax):
    return -alpha*Imax*I/(I+Imax)

# integrate saturation model over a distance d
def saturation(I0, alpha=1, Imax=1, d=np.linspace(0,10,1000)):
    dz=d[1]-d[0]
    I=np.zeros(d.size+1)
    dI=np.zeros(d.size)
    I[0]=I0

    for i in np.arange(1,d.size+1):
        f1 = satModel(I[i-1],alpha,Imax)
        f2 = satModel(I[i-1] + dz * f1 / 2.0,alpha,Imax)
        f3 = satModel(I[i-1] + dz * f2 / 2.0,alpha,Imax)
        f4 = satModel(I[i-1] + dz * f3, alpha,Imax)
        dI[i-1] = (f1 + 2.0 * f2 + 2.0 * f3 + f4 ) / 6.0
        # dI[i-1] = dz * dI[i-1]
        I[i] = I[i-1] + dz * dI[i-1]
    I=I[:-1]

    return I,dI*dz

''''' MODELS FOR ABSORPTION COEFFICIENT '''''
def sqrtDOS(E,Eg=1,alpha0=1):
    loc=np.amin(np.where(E > Eg))
    alpha=np.zeros(E.size)
    alpha[loc:]=np.sqrt(E[loc:]-Eg)
    alpha=alpha/np.amax(alpha)*alpha0

    return alpha
