import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator as plt_aml

def fig(sx=10,sy=7,nMinor=2):
    fig,ax=plt.subplots(figsize=figsize)
    minorTicker(ax,nMinor=nMinor)
    return fig,ax

def fig1(func,sx=10,sy=7,*args):
    fig,ax=plt.subplots(figsize=(sx,sy) )
    minorTicker(ax)
    func(fig,ax,args)

    return fig,ax

def wf1(fig,ax,*args,labX='Time (ps)',labY='Electric Field (kV/cm)'):
    ax.set_xlabel(labX)
    ax.set_ylabel(labY)

    return

def wf2(fig,ax,*args,labX='Time (ps)',labY='Electric Field (a.u)'):
    ax.set_xlabel(labX)
    ax.set_ylabel(labY)

    return

def spec(fig,ax,*args,labX='Frequency (THz)',labY='Amp. (a.u)'):
    ax.set_xlabel(labX)
    ax.set_ylabel(labY)
    ax.set_yscale('log')

    return


def minorTicker(ax,nMinor=2):
    ax.xaxis.set_minor_locator(plt_aml(nMinor))
    ax.yaxis.set_minor_locator(plt_aml(nMinor))

    return 

def colorAxis(ax,cR='k',cL='C1',cT='C0'):
    ax.spines['left'].set_color(cL)
    ax.spines['right'].set_color(cR)
    ax.tick_params(axis='y',which='both',color=cT,labelcolor=cT)
    
    return
