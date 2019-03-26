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
dt = 0.1 #0.01, meshlength of t
npxs = 20 #20, number of mesh points of x
a6i=[60,-2,15,-60,20,30,-3] #a kind of indexs if 6 points scheme
#RK??

def deriv_x(a6i,dx,u):
    '''
v is derivatives u to x by a 6 points Tylor with mesh dx
    '''
    dj=np.arange(-3,3)
    v=np.zeros(len(u))
    for j in dj:
        v+=np.roll(u,dj[j])*a6i[j+4]/(a6i[0]*dx)
        #print(v)
    return v

def RK1(c,dx,dt,t,u):
    '''
u'(x) is given before by deriv_x;
now calcu t with 1-order RK
    '''
    npts=t/dt
    for j in arange(npts):
        u=u-c*deriv_x(a6i,dx,u)*dx/dt # un-up=-c*u'(x)*dx/dt
    return u


#main
if __name__ == '__main__':
    #####################initial:
    dx = 2*np.pi/npxs #npxs=20
    x=np.arange(0,2*np.pi,dx)
    u0=sin(x)

    t1=20 #dt=0.01
    t2=50

    ####################solve at t=20,50, and plt:
    u=RK1(c,dx,dt,t1,u0)
    print(u)

    plt.figure(1)
    plt.plot(x,u0)
    plt.show()

    ###################spectrum analy:
    #v=RK1()#calcu once; or diff1
    v=np.arange(npxs)
    v_=np.fft.ifft(u0) # a N=length(v) cpmplex array? whose k_ is 1,2,...,N
    kr=v_.real*dx
    ki=v_.imag*dx
    alpha=np.arange(npxs)*dx
    print(len(alpha),' ',len(ki))

    plt.figure(3)
    #plt.plot(alpha,ki)
    plt.show()


    '''
    p1 = plt.subplot(2, 1, 1)
    plt.title('errors of frequncy on space') 
    plt.plot(values['x'], alpha, color='b') 
    plt.plot(values['x'], alpha, linewidth=1.5, color='r')
    plt.ylabel('k_i') 
    plt.axis([-0.5, 0.5, -0.05, 1.05])

    ap2 = plt.subplot(2, 1, 2) 
    plt.plot(values['x'], alpha, color='b') 
    plt.plot(values['x'], alpha, linewidth=1.5, color='r')
    plt.ylabel('k_r') 
    plt.xlabel('alpha')
    plt.show()
    '''

#    def invsFT(k):
    '''
inverse Fourier tran;
return the spectrum v_, of a certain freq the wavenumber k, of the difference v
    '''
#    return k;
    #print the values
##### main end
