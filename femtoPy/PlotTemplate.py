from matplotlib.pyplot import subplots as mpl_subplots
from matplotlib.ticker import AutoMinorLocator as plt_aml

def fig(sx=10,sy=7,nMinor=2):
    fig,ax=mpl_subplots(figsize=(sx,sy))
    minorTicker(ax,nMinor=nMinor)
    return fig,ax

def wf1(sx=10,sy=7,nMinor=2,labX='Time (ps)',labY='Electric Field (kV/cm)'):
    fig,ax=fig(sx=sx,sy=sy,nMinor=nMinor)
    ax.set_xlabel(labX)
    ax.set_ylabel(labY)

    return

def wf2(sx=10,sy=7,nMinor=2,labX='Time (ps)',labY='Electric Field (a.u)'):
    fig,ax=fig(sx=sx,sy=sy,nMinor=nMinor)
    ax.set_xlabel(labX)
    ax.set_ylabel(labY)

    return

def spec(sx=10,sy=7,nMinor=2,ylabel='Amp. (a.u)',xlabel='THz'):
    fig,ax=fig(sx=sx,sy=sy,nMinor=nMinor)
    ax.set_ylabel(labY)
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
    return

def PD(sx=10,sy=7,nMinor=2,ylabel='mV',xlabel='ns'):
    fig,ax=fig(sx=sx,sy=sy,nMinor=nMinor)
    ax.set_ylabel('PD Trace ('+ylabel+')')
    ax.set_xlabel('time ('+xlabel+')')
    return fig,ax

def AC(sx=10,sy=7,nMinor=2,ylabel='a.u',xlabel='fs'):
    fig,ax=fig(sx=sx,sy=sy,nMinor=nMinor)
    ax.set_ylabel('AC Trace ('+ylabel+')')
    ax.set_xlabel('time ('+xlabel+')')
    return fig,ax

def minorTicker(ax,nMinor=2):
    ax.xaxis.set_minor_locator(plt_aml(nMinor))
    ax.yaxis.set_minor_locator(plt_aml(nMinor))

    return 

def colorAxis(ax,cR='k',cL='C1',cT='C0'):
    ax.spines['left'].set_color(cL)
    ax.spines['right'].set_color(cR)
    ax.tick_params(axis='y',which='both',color=cT,labelcolor=cT)
    
    return
