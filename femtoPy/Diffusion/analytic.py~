from scipy.special import erfcx

'analytic expression for initial carrier density with monomolecular'
'recombination'
def dist(x,args):
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
