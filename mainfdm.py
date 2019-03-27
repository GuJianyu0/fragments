#### the prog file: only this
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
from scipy.interpolate import spline

'''
1.equation of u=u(x,t):
frac{\partial u}{\partial t}+a*frac{\partial u}{\partial x}=0
boundary condition and meshlength:
x boundary [0,2*pi]; npxs=20; t=20,50; mesh_t=0.01;
set characteristic velocity c=1;
2.use the 5-order upwind scheme to calculate
'''

# the global vars:
c = 1 #characteristic velocity
dt = 0.01 #meshlength of t
npxs = 20 #number of mesh points of x
a6i=[60,-2,15,-60,20,30,-3] #a kind of indexs if 6 points scheme

# called functions
def deriv_x(a6i,dx,u):
    '''
v is derivatives u to x by a 6 points Tylor with mesh dx
    '''
    dj=np.arange(-3,3)
    v=np.zeros(len(u))
    for j in dj:
        v+=np.roll(u,-j)*a6i[j+4]/(a6i[0]*dx) #not roll dj[j], but -j!!!!!!!!!
    return v

def RK3(c,dx,dt,t,u):
    '''
u'(x) is given before by deriv_x;
now calcu t with 1-order RK
    '''
    npts=t/dt
    for j in arange(npts):
        # 1st stage, un-up=-c*u'(x)*dt; 
        u0=u
        dF=c*deriv_x(a6i,dx,u)
        u=u0-dt*dF #if noly RK1, it will diverse!!!!
        # 2nd stage
        dF=c*deriv_x(a6i,dx,u)
        u=0.75*u0+0.25*(u-dt*dF)
        # 3rd stage
        dF=c*deriv_x(a6i,dx,u)
        u=(u0+2*(u-dt*dF))/3
    return u

def exacsolu(c,x,t,u): #!!!!simplifid...
    return sin(x-c*t)

def l2norm(u1,u2):
    return ( sum( (u1-u2)**2 )/len(u1) )**0.5

#main
if __name__ == '__main__':
    ##################### initial, npxs=200, dt=0.01
    dx = 2*np.pi/npxs
    x=np.arange(0,2*np.pi,dx)
    #u01=np.ones(len(x))
    #u01[x>=3]=0
    u0=sin(x)

    t1=20
    t2=50

    #################### solve at t=20,50, and plt:
    u1_exac=exacsolu(c,x,t1,u0)
    u1=RK3(c,dx,dt,t1,u0)
    print('u1 is \n',u1)

    plt.figure(1)
    #plt.plot(x,u1_exac,label='u0',linewidth=2.4)
    x_new = np.linspace(min(x),max(x),50) 
    y_smooth = spline(x, u1_exac, x_new) #should exac solu be smoothed ?????
    plt.plot(x_new, y_smooth)
    plt.plot(x,u1,label='u1',color='r',linewidth=0.8) # a less wider line to see
    plt.legend()
    plt.xlabel('x') 
    plt.ylabel('u(t=20)')
    
    u2_exac=exacsolu(c,x,t2,u0)
    u2=RK3(c,dx,dt,t2,u0)
    plt.figure(2)
    #plt.plot(x,u2_exac,label='u0',linewidth=2.4)
    x_new = np.linspace(min(x),max(x),50) 
    y_smooth = spline(x, u2_exac, x_new) #should exac solu be smoothed ?????
    plt.plot(x_new, y_smooth)
    plt.plot(x,u2,label='u2',color='r',linewidth=0.8)
    plt.legend()
    plt.xlabel('x') 
    plt.ylabel('u(t=50)')

    ################### spectrum analy:
    X=np.arange(npxs)
    i_a=np.arange(-3,3)
    alpha=X*dx
    v=deriv_x(a6i,dx,u0)

    #vcap, spectrum of v
    vcap=np.zeros(npxs)+0j
    #vcap_ifft=np.fft.ifft(v) #fft????????
    for k in X:
        #for j in X:  #?????????????????????????????????
        #    vcap[k]+=v[j]*exp(-1j*k*x[j])/npxs
        for i in i_a:
            vcap[k]+=exp(1j*k*dx*i)*a6i[i+4]/a6i[0]/dx
    print('vcap is \n',vcap)

    #kcap=kr+ki*1j, and plot
    kr=vcap.real*dx
    ki=vcap.imag*dx

    npxs2=npxs//2
    plt.figure(3)
    plt.plot(alpha[:npxs2],alpha[:npxs2])
    plt.plot(alpha[:npxs2],ki[:npxs2],label='ki')
    plt.legend()
    plt.xlabel('alpha') 
    plt.ylabel('ki')

    plt.figure(4)
    plt.plot(alpha[:npxs2],alpha[:npxs2])
    plt.plot(alpha[:npxs2],kr[:npxs2],label='kr')
    plt.legend()
    plt.xlabel('alpha') 
    plt.ylabel('kr')

    plt.show()

    ################# print l2-mod error
    u1_l2error=l2norm(u1,u1_exac)
    u2_l2error=l2norm(u2,u2_exac)
    print(u1_l2error,u2_l2error)
##### end main
