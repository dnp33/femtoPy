from matplotlib.pyplot import subplots as mpl_subplots
from matplotlib.ticker import AutoMinorLocator as plt_aml
from numpy import array as np_array

def figure(sx=10,sy=7,nMinor=1):
    fig,ax=mpl_subplots(figsize=(sx,sy))
    minorTicker(ax,nMinor=nMinor)
    return fig,ax

def wf1(sx=10,sy=7,nMinor=1,xlabel='Time (ps)',ylabel='Electric Field (kV/cm)'):
    fig,ax=figure(sx=sx,sy=sy,nMinor=nMinor)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    return fig,ax

def wf2(sx=10,sy=7,nMinor=1,xlabel='Time (ps)',ylabel='Electric Field (a.u)'):
    fig,ax=figure(sx=sx,sy=sy,nMinor=nMinor)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    return fig,ax

def spec(xlabel='THz',sx=10,sy=7,nMinor=1,ylabel='Amp. (a.u)'):
    fig,ax=figure(sx=sx,sy=sy,nMinor=nMinor)
    ax.set_ylabel(ylabel)
    if xlabel=='THz':
        ax.set_xlabel('Frequency (THz)')
    elif xlabel=='eV':
        ax.set_xlabel('Energy (eV)')
    elif xlabel=='meV':
        ax.set_xlabel('Energy (meV)')
    elif xlabel=='nm':
        ax.set_xlabel('Wavelength (nm)')
    elif xlabel=='um':
        ax.set_xlabel(r'Wavelength ($\mu$m)')
    else:
        ax.set_xlabel(xlabel)
    return fig,ax

def PD(sx=10,sy=7,nMinor=1,ylabel='mV',xlabel='ns'):
    fig,ax=figure(sx=sx,sy=sy,nMinor=nMinor)
    ax.set_ylabel('PD Trace ('+ylabel+')')
    ax.set_xlabel('time ('+xlabel+')')
    return fig,ax

def AC(sx=10,sy=7,nMinor=1,ylabel='a.u',xlabel='fs'):
    fig,ax=figure(sx=sx,sy=sy,nMinor=nMinor)
    ax.set_ylabel('AC Trace ('+ylabel+')')
    ax.set_xlabel('time ('+xlabel+')')
    return fig,ax

def minorTicker(ax,nMinor=1):
    axType=type(ax)
    if axType==type(np_array([])) or axType==type([]):
        for i in range(len(ax)):
            ax[i].xaxis.set_minor_locator(plt_aml(nMinor+1))
            ax[i].yaxis.set_minor_locator(plt_aml(nMinor+1))
    else:
        ax.xaxis.set_minor_locator(plt_aml(nMinor+1))
        ax.yaxis.set_minor_locator(plt_aml(nMinor+1))

    return 

def colorAxis(ax,cR='k',cL='C1',cT='C0'):
    ax.spines['left'].set_color(cL)
    ax.spines['right'].set_color(cR)
    ax.tick_params(axis='y',which='both',color=cT,labelcolor=cT)
    
    return
