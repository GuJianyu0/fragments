#Reference: the program referenced ""
#### the main funtion file
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
#from matplotlib.ticker import MultipleLocator, FormatStrFormatter
#import 

'''equation: 
frac{\partial u}{\partial t}+a*frac{\partial u}{\partial x}=0
boundary condition:
x boundary [0,2*pi]; npxs=20; t=20,50; mesh_t=0.01;
set characteristic velocity c=1;
'''

# the global vars:
c = 1 #characteristic velocity
dt = 0.01 #0.01, meshlength of t
npxs = 200 #20, number of mesh points of x
a6i=[60,-2,15,-60,20,30,-3] #a kind of indexs if 6 points scheme

# called functions
def deriv_x(a6i,dx,u):
    '''
v is derivatives u to x by a 6 points Tylor with mesh dx
    '''
    dj=np.arange(-3,3)
    v=np.zeros(len(u))
    for j in dj:  #not roll dj[j], but -j!!!!!!!!!
        v+=np.roll(u,-j)*a6i[j+4]/(a6i[0]*dx) #v=(-2*np.roll(u,3)+15*np.roll(u,2)-60*np.roll(u,1)+20*np.roll(u,0)+30*np.roll(u,-1)-3*np.roll(u,-2))/(60*dx)
    return v

def RK(c,dx,dt,t,u):
    '''
u'(x) is given before by deriv_x;
now calcu t with 1-order RK
    '''
    npts=t/dt
    for j in arange(npts):
        # 1st stage  # un-up=-c*u'(x)*dt; if RK1, diverse!!!!
        u0=u
        dF=c*deriv_x(a6i,dx,u)
        u=u0-dt*dF
        # 2nd stage
        dF=c*deriv_x(a6i,dx,u)
        u=0.75*u0+0.25*(u-dt*dF)
        # 3rd stage
        dF=c*deriv_x(a6i,dx,u)
        u=(u0+2*(u-dt*dF))/3
    return u

#main prog
if __name__ == '__main__':
    #####################initial:
    dx = 2*np.pi/npxs #npxs
    x=np.arange(0,2*np.pi,dx)
    u01=np.ones(len(x))
    u01[x>=3]=0
    u0=sin(x)

    t1=20 #dt=0.01
    t2=50

    ####################solve at t=20,50, and plt:
    u1=RK(c,dx,dt,t1,u0)
    #print(u0,'\n',u1)

    plt.figure(1)
    plt.plot(x,u1)
    plt.show()

    ###################spectrum analy:
    X=np.arange(npxs)
    alpha=X*dx
    v=deriv_x(a6i,dx,u0)

    #vcap, spectrum of v
    vcap_ifft=np.fft.ifft(v) #? a N=length(v) cpmplex array, whose k_ is 1,2,...,N
    vcap=np.zeros(npxs)+0j
    for k in X:
        for j in X:
            vcap[k]=v[j]*exp(-1j*k*x[j])/npxs
    print(vcap)

    #kcap=kr+ki*1j
    kr=vcap.real/dx #_ifft
    ki=vcap.imag/dx

    plt.figure(3)
    plt.plot(x,x)
    plt.plot(alpha,ki)#/15e-5
    plt.show()

    #print the values
##### end
