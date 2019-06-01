import matplotlib.gridspec as gridspec
import matplotlib.ticker as FEMTOTICKER
import matplotlib.pyplot as plt

import numpy as np
import numpy.fft as fft

import scipy.optimize as opt

'make a figure with the amplitude and phase for the THz pulse'
def fftFig(figsize=(10,8),left=.13,right=.87,mid=.7,top=.9,bottom=.13,hspace=0):
    fig=plt.figure(figsize=figsize)
    'make amplitude figure'
    gs2=gridspec.GridSpec(1,1)
    gs2.update(left=left,right=right,top=mid,bottom=bottom,hspace=0)
    axA=plt.subplot(gs2[0,0])
    axA.set_xlabel('frequency (THz)')
    axA.set_ylabel('amplitude (a.u)')
    axA.xaxis.set_minor_locator(FEMTOTICKER.AutoMinorLocator(2))
    axA.yaxis.set_minor_locator(FEMTOTICKER.AutoMinorLocator(2))

    'make phase figure'
    gs1=gridspec.GridSpec(1,1)
    gs1.update(left=left,right=right,top=top,bottom=mid)
    axP=plt.subplot(gs1[0,0])
    axP.set_ylabel('phase\n(a.u)')
    # axP.set_xticks([])
    # axP.set_yticks([])
    axP.xaxis.set_minor_locator(FEMTOTICKER.AutoMinorLocator(2))
    axP.yaxis.set_minor_locator(FEMTOTICKER.AutoMinorLocator(2))

    return fig,axA,axP

##################################
#This function will determine the
#thickness of a sample given
#a range of thicknesses assumed
#in addition to the best fit
#index and absorption coefficent
#according to the parameter
#extraction model from Koch 07a OE
#
#INPUT:
#   rdat -> reference data (1D array?)
#   sdat -> sample data (1D array?)
#   tpts -> time data (seconds, 1D array?)
#   fmin -> minimum frequency for evaluating n, a and t (Hz)
#   fmax -> max frequency
#   tmin -> minimum thickness to consider (m)
#   tmax -> maximum thickness to consider (m)
#   dt -> thickness step size to evaluate
#OUTPUT:
#   nb -> real part of the index of refraction from best fit to thickness
#   kb -> imaginary part of the index of refraction from best fit to
#   thickness
#   tb -> best fit thickness
#   TV -> total variation measure (which should have nice minimum if the
#   extraction algorithm worked).
##################################




#compute complex index and thickness following Matt Reid from UNBC, adapted from Koch 07a OE
def compute_n_a_t(rdat,sdat,tpts,fmin,fmax,tmin,tmax,dt):
    #transfer function
    def H_theory(nc,fp,thick,nrefl):
        c0 = 2.99796e8 #speed of light
        n0 = 1.00027   #index of air
        T01 = 2*n0/(n0+nc) #transmission coefficient?
        T10 = 2*nc/(nc+n0) #transmission coefficient?
        R10 = (n0-nc)/(n0+nc) #reflection coefficient?
        P0 = np.exp(1j*(2*np.pi*fp)*n0*thick/c0) #??
        #%size(fp)
        #%size(nc)
        P1 = np.exp(-1j*(2*np.pi*fp)*nc*thick/c0) #??
        
        sumdelta = 0  #??
        for ii in np.arange(0,nrefl):
            sumdelta = sumdelta+(R10**2 * P1**2)**ii

            HOut = P0*T01*P1*T10*(1+sumdelta)
            
        return HOut
    #constants
    c = 2.99796e8 #speed of light
    nair = 1.00027 #index of air
    thickness = 1./2.*(tmin+tmax) # initial guess at thickness

    #remove offset
    rdat = rdat-np.average(rdat)   #remove dc component of reference waveform
    sdat = sdat - np.average(sdat) #remove dc component of sample waveform

    #compute FFT's
    fr = fft.fftshift(fft.fft(rdat))
    fs = fft.fftshift(fft.fft(sdat))
    fpr = fft.fftshift(fft.fftfreq(tpts.size,tpts[2]-tpts[1]))

    #locate fmax and min and trim data
    lmax = np.amin(np.where((fpr-fmax)>0))
    lmin = np.amax(np.where((fpr-fmin)<0))  
    fpr = fpr[lmin:lmax]                    
    fr = fr[lmin:lmax]                      
    fs = fs[lmin:lmax]                      

    #guess for n and k using the time delay total decrease in amp
    mr = np.amax(rdat)        # max value of EO signal
    lr = np.where(rdat==mr)   # index of max
    trm = tpts[lr]            # time of max
    
    ms = np.amax(sdat)      
    lss = np.where(sdat==ms)
    tsm = tpts[lss]
    
    delta_t = tsm - trm

    
    n_approx = np.zeros(fpr.size)+1
    n_approx = (delta_t*c/thickness + nair)*n_approx 

    # guess k based on absorption
    k_approx =-1/thickness*c/(2*np.pi*fpr)*np.log(np.absolute(fs)/np.absolute(fr)) 

    ####################################################
    #calculate numer of multiple reflections to include#
    ####################################################
    t_refl = 2*thickness*n_approx[1]/c            # time between reflections???
    n_refl = np.floor((tpts[tpts.size-1]-tpts[1])/t_refl)    # number of reflections to consider based on initial guess?

    ##############################
    #Compute thickness
    #by minimizing total
    #variation between ref. and trans.
    ##############################

    N=np.round((tmax-tmin)/tstep+1)
    thpts = np.linspace(tmin,tmax,N)

    #opts=optimset('MaxFunEvals',1e7,'MaxIter',1e6)   #parameters for fminsearch
    H_exp = fs/fr   #experimental transfer function
    
    # functions defining the difference between theoretical transfer function, and measured transfer function
    # I have no idea where the theoretical transfer function comes from, this is a blackbox to me...
    def Mw(npred,kpred):
        return np.absolute(H_theory(npred-1j*kpred,fpr,thickness,n_refl)) - np.absolute(H_exp)
    def Aw(npred,kpred):
        return np.unwrap(np.angle(H_theory(npred-1j*kpred,fpr,thickness,n_refl))) - np.unwrap(np.angle(H_exp))
    def Err2Min(x0):
        x0=np.reshape(x0,(x0.size/2,2))
        return np.sum(np.absolute(Mw(x0[:,0],x0[:,1]))+np.absolute(Aw(x0[:,0],x0[:,1])))
    

    
    # will hold optimized values of n and k for each thickness guess
    ncalc = np.empty((n_approx.size,thpts.size))
    kcalc = np.empty((n_approx.size,thpts.size))
    
    # this is an estimate of the error from H_exp to H_theory
    TV = np.empty(thpts.size)
    
    ##################################################
    # loop over each thickness and calculate optimum #
    # parameters of n and k which minimize difference#
    # between experimental and theoretical transfer  #
    # function                                       #
    ##################################################
    for ii in np.arange(0,thpts.size):
        print('step ',ii+1,' out of ',thpts.size)
        thickness = thpts[ii]          # thickness of current iteration
        x0=np.empty((n_approx.size,2))   # 2d array with n and k

        x0[:,0]=n_approx               
        x0[:,1]=k_approx
        
        x0=np.ndarray.flatten(x0)
        
        ##########################################
        # minimize difference between theoretical#
        # and experimental transfer function as  #
        # a function of n(f), k(f): this is where#
        # the dirty work is done!                #
        ##########################################
        Soln = opt.fmin(Err2Min,x0,maxfun=1e7,maxiter=1e6)
        Soln = np.reshape(Soln,(Soln.size/2,2))
        
        #update guesses for n and k
        ncalc[:,ii] = Soln[:,0]
        kcalc[:,ii] = Soln[:,1]
        
        # the fuck is going on here?
        Dm = np.absolute(np.roll(ncalc[:,ii],1)-ncalc[:,ii])+np.absolute(np.roll(kcalc[:,ii],1)-kcalc[:,ii])
        TV[ii] = np.sum(Dm[0:Dm.size])

    tvmin = np.amin(TV)
    ltvmin = np.where(TV==tvmin)
    tb = thpts[ltvmin]
    
    nb = ncalc[:,ltvmin]
    nb = np.reshape(nb,nb.size)
    
    kb = kcalc[:,ltvmin]
    kb = np.reshape(kb,kb.size)

    return nb, kb, tb, TV, fpr
