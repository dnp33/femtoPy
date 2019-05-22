import matplotlib.pyplot as plt

def fig(sx=10,sy=7):
    fig,ax=plt.subplots(figsize=figsize)
    ax.xaxis.set_minor_locator(FEMTOTICKER.AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(FEMTOTICKER.AutoMinorLocator(2))

    return fig,ax

def fig1(func,sx=10,sy=7,*args):
    fig,ax=plt.subplots(figsize=(sx,sy))
    ax.xaxis.set_minor_locator(FEMTOTICKER.AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(FEMTOTICKER.AutoMinorLocator(2))

    fig,ax=func(fig,ax,*args)

    return fig,ax

def wf1(fig,ax,labX='Time (ps)',labY='Electric Field (kV/cm)'):
    ax.set_xlabel(labX)
    ax.set_ylabel(labY)

    return fig,ax

def wf2(fig,ax,labX='Time (ps)',labY='Electric Field (a.u)'):
    ax.set_xlabel(labX)
    ax.set_ylabel(labY)

    return fig,ax

def spec(fig,ax,labX='Frequency (THz)',labY='Amp. (a.u)'):
    ax.set_xlabel(labX)
    ax.set_ylabel(labY)
    ax.set_yscale('log')

    return fig,ax
