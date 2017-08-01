import numpy as np  

'class to hold the grid parameters for the FDTD code'
class gridParams:
    def __init__(self,dt=0.05,dy=0.01,y_min=0,y_max=5,t_min=0,t_max=100):
        'grid parameters'
        self.dt = dt # grid size for time
        self.dy = dy # grid size for space
        self.y = np.arange(y_min,y_max+dy,dy) # y-array for grid
        self.t = np.arange(t_min,t_max+dt,dt) # t-array for grid
        
        return

