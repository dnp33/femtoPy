import numpy as np
import femtoPy.PlotTemplate as ptemp
import matplotlib as matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as FEMTOTICKER

np.seterr(all='raise')

plt.rcParams['figure.figsize']=[10,7]

plt.rcParams['axes.labelsize']=23
plt.rcParams['axes.titlesize']=30
#plt.rcParams['axes.grid']=True
plt.rcParams['axes.linewidth']=3
plt.rcParams['axes.xmargin']=0
plt.rcParams['axes.ymargin']=.01

plt.rcParams['lines.linewidth']=3

plt.rcParams['xtick.labelsize']=23
plt.rcParams['ytick.labelsize']=23

plt.rcParams['xtick.major.size']=12
plt.rcParams['xtick.minor.size']=6
plt.rcParams['xtick.major.width']=2
plt.rcParams['xtick.minor.width']=2
plt.rcParams['xtick.minor.visible']=True

plt.rcParams['ytick.major.size']=12
plt.rcParams['ytick.minor.size']=6
plt.rcParams['ytick.major.width']=2
plt.rcParams['ytick.minor.width']=2
plt.rcParams['ytick.minor.visible']=True

plt.rcParams['xtick.direction']='out'
plt.rcParams['ytick.direction']='out'

plt.rcParams['legend.fontsize']=20
plt.rcParams['legend.fancybox']=True
plt.rcParams['legend.scatterpoints']=1

def minorTicker(ax):
    ax.xaxis.set_minor_locator(FEMTOTICKER.AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(FEMTOTICKER.AutoMinorLocator(2))

def colorAxis(ax,cR='k',cL='C1',cT='C0'):
    ax.spines['left'].set_color(cL)
    ax.spines['right'].set_color(cR)
    ax.tick_params(axis='y',which='both',color=cT,labelcolor=cT)
    
    return

"""
from imp import find_module
from imp import load_module

from os.path import expanduser
home=expanduser("~")

fp,pathname,description=find_module('timeResolved',[home+'/Documents/python'])
TR=load_module('timeResolved',fp,pathname,description)
"""
