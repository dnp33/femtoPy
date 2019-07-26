"""
this module contains a set of default plots for commonly used formats

Functions
---------
figure : makes a 10x7 figure
fig2col : a 20x7 two column figure

twinx_ax : twinx & color/label an axis
minorTicker : sets the number of minor ticks
colorAxis : colors the spines, ticks, and labels of an axis
colorAxis2 : " " for a twinned figure

sigma : conductivity figure. optionally twins
epsilon : dielectric function figure. optionally twins
index : index figure. optionally twins
wf1/wf2 : THz waveform figures
spec : spectrum figure
voltage : voltage figure
AC : autocorrelation trace figure.
"""
from matplotlib.pyplot import subplots as mpl_subplots
from matplotlib.ticker import AutoMinorLocator as plt_aml
from numpy import array as np_array

def figure(sx=10,sy=7,nMinor=1,fmt=None,**kwargs):
    """
    Returns a figure of size sx by sy (default 10x7) with a single axis.

    Parameters
    ----------
    sx (sy) : width (height)
    nMinor : number of minor ticks (default 1)
    fmt : formatting function (optional)
    **kwargs : refer to fmt function
    """
    fig,ax=mpl_subplots(figsize=(sx,sy))
    minorTicker(ax,nMinor=nMinor)
    if type(fmt) != type(None):fig,ax=fmt(fig=fig,ax=ax,**kwargs)
    return fig,ax

def fig2col(sx=20,sy=7,nMinor=1,fmt=None,**kwargs):
    """
    returns a figure of size sx by sy with 2 axes
    
    Parameters
    ----------
    sx (sy) : width (height)
    nMinor : number of minor ticks
    fmt : formatting function (optional)
    **kwargs : see fmt function
    """
    fig,(ax0,ax1)=mpl_subplots(figsize=(sx,sy),ncols=2)
    minorTicker([ax0,ax1])
    if type(fmt) !=type(None):ax0,ax1=fmt(fig,ax,**kwargs)
    return fig,ax0,ax1

def twinx_ax(ax,nMinor=1,cL='C0',cR='C1',xlabel='',Rlabel='',Llabel=''):
    """
    put a second y axis on an axis

    Parameters
    ----------
    ax : axis to be twinned
    nMinor : number of minor ticks for twin axis (default 1)
    cL(cR) : color of left(right) axis spine, ticks, and label, default C0(C1)
    xlabel : x axis label (default '')
    Rlabel (Llabel) : right/left y axis labels
    """
    ax0=ax.twinx(); minorTicker(ax0,nMinor=nMinor)
    ax,ax0=colorAxis2(ax,ax0,cR=cR,cL=cL)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(Llabel)
    ax0.set_ylabel(Rlabel)
    return ax,ax0

def sigma(fig=None,ax=None,sx=10,sy=7,nMinor=1,xlabel='frequency (THz)',
          twin=True,units='$S$ $cm^{-1}$',cR='C1',cL='C0'):
    """
    conductivity vs frequency plot

    Parameters
    ----------
    fig,ax : default None (makes new figure)
    sx(sy) : width (height)
    nMinor : number of minor ticks
    xlabel : default - frequency (THz)
    units : default - S cm^-1
    cR(cL) : color of right(left) axis/ticks/labels
    nMinor : default 1
    """
    if type(fig)==type(None):
        fig,ax=figure(sx=sx,sy=sy,nMinor=nMinor)
    if twin:
        return fig,twinx_ax(ax,cR=cR,cL=cL,xlabel=xlabel,
                       Rlabel=r'Im{$\sigma$} ('+units+')',
                       Llabel=r'Re{$\sigma$} ('+units+')')
    else:
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r'$\tilde{\sigma}}$ ('+units+')')
        return fig,ax

def epsilon(fig=None,ax=None,sx=10,sy=7,nMinor=1,xlabel='frequency (THz)',
            twin=True,units='$F$ $m^{-1}$',cR='C1',cL='C0'):
    """
    conductivity vs frequency plot

    Parameters
    ----------
    fig,ax : default None (makes new figure)
    sx(sy) : width (height)
    nMinor : number of minor ticks
    xlabel : default - frequency (THz)
    units : default - F m^-1
    cR(cL) : color of right(left) axis/ticks/labels
    nMinor : default 1
    """
    if type(fig)==type(None):
        fig,ax=figure(sx=sx,sy=sy,nMinor=nMinor)
    if twin:
        return fig,twinx_ax(ax,cR=cR,cL=cL,xlabel=xlabel,
                       Rlabel=r'Im{$\epsilon$} ('+units+')',
                       Llabel=r'Re{$\epsilon$} ('+units+')')
    else:
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r'$\tilde{\epsilon}}$ ('+units+')')
        return fig,ax

def index(fig=None,ax=None,sx=10,sy=7,nMinor=1,xlabel='frequency (THz)',
          twin=True,cR='C1',cL='C0',y2=r'$\alpha (cm^{-1})$'):
    """
    conductivity vs frequency plot

    Parameters
    ----------
    fig,ax : default None (makes new figure)
    sx(sy) : width (height)
    nMinor : number of minor ticks
    xlabel : default - frequency (THz)
    which : extinction coefficient (k) or absorption coefficient (alpha)
    cR(cL) : color of right(left) axis/ticks/labels
    nMinor : default 1
    """
    if type(fig)==type(None):
        fig,ax=figure(sx=sx,sy=sy,nMinor=nMinor)
    if twin:
        return fig,twinx_ax(ax,cR=cR,cL=cL,xlabel=xlabel,
                        Llabel='n',Rlabel=y2)
    else:
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r'$\tilde{n}$')
        return fig,ax

def wf1(fig=None,ax=None,sx=10,sy=7,nMinor=1,xlabel='time (ps)',
        ylabel='electric field (kV/cm)'):
    """
    waveform plot

    Parameters
    ---------
    fig,ax : default None (makes new figure)
    sx(sy) : width (height)
    nMinor : number of minor ticks
    xlabel, ylabel : default time (ps) & E-field (kV/cm)
    """
    if type(fig)==type(None):
        fig,ax=figure(sx=sx,sy=sy,nMinor=nMinor)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    return fig,ax

def wf2(fig=None,ax=None,sx=10,sy=7,nMinor=1,xlabel='time (ps)',
        ylabel='E-field (a.u)'):
    """
    waveform plot

    Parameters
    ----------
    fig,ax : default None (makes new figure)
    sx(sy) : width (height)
    nMinor : number of minor ticks
    xlabel, ylabel : default time (ps) & E-field (a.u)
    """
    if type(fig)==type(None):
        fig,ax=figure(sx=sx,sy=sy,nMinor=nMinor)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    return fig,ax

def spec(fig=None,ax=None,sx=10,sy=7,nMinor=1,label='frequency',
         units='THz',ylabel='mp. (a.u)'):
    """
    spectrum plot
    
    Parameters
    ----------
    fig,ax : default None (makes new figure)
    sx(sy) : width (height)
    nMinor : number of minor ticks
    label : x axis parameter (default - frequency)
    units : units of x axis (default THz)
    ylabel : default amp. (a.u)
    """
    if type(fig)==type(None):
        fig,ax=figure(sx=sx,sy=sy,nMinor=nMinor)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel+' ('+xunits+')')

    return fig,ax

def voltage(fig=None,ax=None,sx=10,sy=7,nMinor=1,yunits='mV',xunits='ns'):
    """
    voltage vs time plot

    Parameters
    ----------
    fig,ax : default None (makes new figure)
    sx(sy) : width (height)
    nMinor : number of minor ticks
    yunits : default mV
    xunits : default ns
    """
    if type(fig)==type(None):
        fig,ax=figure(sx=sx,sy=sy,nMinor=nMinor)
    ax.set_ylabel('voltage ('+ylabel+')')
    ax.set_xlabel('time ('+xlabel+')')
    return fig,ax

def AC(fig=None,ax=None,sx=10,sy=7,nMinor=1,ylabel='AC trace (a.u)',
       xlabel='time (fs)'):
    """
    autocorrelation vs time graph

    Parameters
    ----------
    fig,ax : default to None (makes a figure)
    xlabel (ylabel) : default time (fs) & AC trace (a.u)
    """
    if type(fig)==type(None):
        fig,ax=figure(sx=sx,sy=sy,nMinor=nMinor)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    return fig,ax

def minorTicker(ax,nMinor=1):
    """sets the number of minor ticks to nMinor"""
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
    """
    colors the y axis

    Parameters
    ----------
    cR : right axis spine color
    cL : left axis spine color
    cT : tick color 
    """
    
    ax.spines['left'].set_color(cL)
    ax.spines['right'].set_color(cR)
    ax.tick_params(axis='y',which='both',color=cT,labelcolor=cT)
    ax.yaxis.label.set_color(cT)
    
    return ax

def colorAxis2(ax0,ax1,cR='C0',cL='C1'):
    """
    same as color axis, meant for a twinned figure
    """
    ax0.spines['left'].set_color(cL)
    ax0.spines['right'].set_color(cR)
    ax0.tick_params(axis='y',which='both',color=cL,labelcolor=cL)
    ax0.yaxis.label.set_color(cL)
    
    ax1.spines['left'].set_color(cL)
    ax1.spines['right'].set_color(cR)
    ax1.tick_params(axis='y',which='both',color=cR,labelcolor=cR)
    ax1.yaxis.label.set_color(cR)

    return ax0,ax1
