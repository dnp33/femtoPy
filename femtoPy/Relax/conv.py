"""
code to efficiently calculate the convolution of a response function and a 
fit function
"""

from numpy import exp as np_exp, sqrt as np_sqrt
from numpy.fft import rfft as fft_rfft, irfft as fft_irfft

def gaussResponse(t,sigma):
    return np_exp(-(t/sigma)**2/2)/(sigma*np.sqrt(2*np.pi))

def convCalc(t,dat,resp):
    dat=fft_rfft(dat); resp=fft_rfft(resp)
    dat=fft_irfft(dat*resp)

    return dat*(t[1]-t[0])
