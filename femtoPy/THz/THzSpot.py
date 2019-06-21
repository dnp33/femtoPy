# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 14:08:53 2019

@author: cam_hough21
"""

from numpy import loadtxt as np_loadtxt;
import matplotlib.pyplot as plt

class THzSpot:
    def __init__(self, filename1='spot',filename2='bg'):
        self.im=np_loadtxt(filename1)
        self.bg=np_loadtxt(filename2)
        self.spot=self.im-self.bg
        
        return
    
    def plot(self):
        fig,ax=plt.subplots()
        ax.contourf(self.spot)
        
        return fig,ax
    
