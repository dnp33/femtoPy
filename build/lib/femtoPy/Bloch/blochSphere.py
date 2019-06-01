from femtoPy.imports import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as anim

def f(R):
    k=np.cross(R,Omega)
    k[0]=k[0]-Gamma*R[0]/2
    k[1]=k[1]-Gamma*R[1]/2
    k[2]=k[2]+Gamma*(1-R[2])

    return k

def BlochEq(R,Omega,dt):
    k1=f(R)*dt
    k2=f(R+k1/2)*dt
    k3=f(R+k2/2)*dt
    k4=f(R+k3)*dt
    
    R=R+(k1+2*(k2+k3)+k4)/6

    return R

R=np.array([0,0,1])
dt=0.01
t=np.zeros(1)

n=np.zeros(1)

RabiFreq=8
delta=0.
Omega=np.array([RabiFreq,0,delta])

Gamma=0.1

trace=np.array([0,0,1])


#animation
def drawSingle(ax):
    global R,trace,t,Omega,n
    t=np.append(t,t[-1]+dt)
    if t[-1] > np.pi/RabiFreq:
        Omega[0]=0
        
    R=BlochEq(R,Omega,dt)
    trace=np.append(trace,R)
    n=np.append(n,(1-R[2])  **2)
    ax.plot_wireframe(x,y,z,color='k',alpha=0.2,linewidth=0.3)
    ax.plot([0,R[0]],[0,R[1]],[0,R[2]],color='b')

    line=np.reshape(trace,(int(trace.size/3),3))
    ax.plot(line[:,0],line[:,1],line[:,2],linewidth=0.5,color='k')
    
    ax.set_xlim(-Max,Max)
    ax.set_ylim(-Max,Max)
    ax.set_zlim(-Max,Max)

    # get rid of panes
    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

    # Get rid of the spines
    ax.w_xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
    ax.w_yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
    ax.w_zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))

    # Get rid of the ticks                          
    ax.set_xticks([])                               
    ax.set_yticks([])                               
    ax.set_zticks([])

    ax.text(1,1,1,'|R|='+str(np.sqrt(np.sum(R*R))))
    ax.text(1,1,0.89,'t='+str(t[-1]))
    
    return

#set up plots
fig=plt.figure(figsize=(10,10))
ax=fig.add_subplot(111,projection='3d')

u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x = np.cos(u)*np.sin(v)
y = np.sin(u)*np.sin(v)
z = np.cos(v)

Max=1.2
drawSingle(ax)
# ax.plot_surface(x,y,z,color='k',alpha=0.05)
# ax.plot([0,0],[0,0],[0,1])



def animate(phi):
    ax.cla()
    drawSingle(ax)

ani = anim.FuncAnimation(fig,animate,interval=0,blit=False)

plt.show()

fig,ax=figure()
ax.plot(t,n)
ax.set_xlabel('time')
ax.set_ylabel('|c2| squared')

ax.set_ylim(-0.1,1.1)

plt.show()
