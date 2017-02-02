import numpy as np
import matplotlib as matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker

def setTicks(ax):
    locator=matplotlib.ticker.AutoMinorLocator(n=2)
    ax.xaxis.set_minor_locator(locator)
    ax.yaxis.set_minor_locator(locator)

def newFig(x=10,y=8):
    fig,ax=plt.subplots(figsize=(x,y))
    setTicks(ax)

    return fig,ax

plt.rcParams['figure.figsize']=(10,8)

plt.rcParams['axes.labelsize']=23
plt.rcParams['axes.titlesize']=30
plt.rcParams['axes.grid']=True
plt.rcParams['axes.linewidth']=3

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


#from imp import find_module
#from imp import load_module

#from os.path import expanduser
#home=expanduser("~")

#fp,pathname,description=find_module('timeResolved',[home+'/Documents/python'])
#TR=load_module('timeResolved',fp,pathname,description)
