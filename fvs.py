# the main prog file
import exact
import numpy as np
import matplotlib.pyplot as plt

#################### instructions
'''
1.Equation of U=U(x,t), U=[rho, u, p]^T:
frac{\partial U}{\partial t}+*frac{\partial f(U)}{\partial x}=0

2.Boundary condition and meshlength:
x is domain[-0.5,0.5]; npxs=200;
t1=0.14; mesh_t=0.001;

3.Use the FVS(Flow Vector Splite) scheme to calculate:
flow: f(U) = A U = S^{-1} \Lamda S U,
where \Lamda is eigenvalue matrix of A, \lamda_1=u, \lamda_2=u-c, \lamda3=u+c;
Steager-Warming scheme splite;
then RK the time direction;

4.Notice each var is a vector in domain of x in order to calculate together (,so that infact U.shape is (3,len(x)).
'''

################### global vars
gamma = 1.4 #diabatic exponent
dt = 0.001 #meshlength of t
npxs = 200 #number of mesh points of x
epsilon=1.e-6 #to avoid 1/0 in lamda+- ???

################### called functions and intervars
def lamda(U,c):
    '''
    eigenvalues of A;
    notice here l[0] represents \lamda_1; lamda_p represent \lamda^{+};
    to return lamda_p, lamda_m
    '''
    l=np.zeros(U.shape)
    u=U[1]
    l[0]=u; l[1]=u-c; l[2]=u+c
    #l=abs(u)+c #LF splite
    lamda_p=(l+np.sqrt(l**2+epsilon**2))/2
    lamda_m=(l-np.sqrt(l**2+epsilon**2))/2
    return lamda_p, lamda_m

def c_sound(U, gamma): #sound speed
    rho=U[0]; p=U[2];
    #p[p<0]=0 #?????
    return np.sqrt(gamma * p/ rho)

def wf(lamda,c,gamma): # just an intervar of flow f 
    return (3.-gamma)*(lamda[1]+lamda[2])*c**2/2/(gamma-1)

def flow(U,lamda,c,gamma):
    f=np.zeros(U.shape)
    rho=U[0]; u=U[1];
    w=wf(lamda,c,gamma)
    f[0]=2*(gamma-1.)*lamda[0]+lamda[1]+lamda[2]
    f[1]=2*(gamma-1.)*lamda[0]*u+lamda[1]*(u-c)+lamda[2]*(u+c)
    f[2]=(gamma-1.)*lamda[0]*u**2+lamda[1]/2*(u-c)**2+lamda[2]/2*(u+c)**2+w
    return f*rho/2/gamma

def bondary_roll2c(a,n,axis=1):#.........
    b=np.roll(a,n,1)
    if n>=0: #translate right n
        for i in np.arange(n):
            b[:,i]=a[:,0]
    if n<0: #translate left |n|
        for i in np.arange(n,0):
            b[:,i]=a[:,-1]
    return b

def deriv_x(f,pm,dx):#bondary, change roll f(U)'s comeflow??????
    '''
    v(i.e. fp_x) is partial derivatives of f to dx, a 6 points upwind scheme with p(plus corresponed to f(lamda+)) to a(characteristic velocity)>0 and m(-) to a<0
    '''
    v=np.zeros(f.shape)
    a6i=[60,-2,15,-60,20,30,-3]
    J=np.arange(-3,3)
    if pm==1: # lamda+, a>0 upwind scheme
        for j in J:
            v+=bondary_roll2c(f,-j,1)*a6i[j+4]/(a6i[0]*dx)
        #v=-(bondary_roll2c(f,1,1)-f)/dx # test: a same error fig as the 5-order

    else: #a<0; not safe, there is no lamda+- but f+-, extra input needed......
        for j in J:
            v-=bondary_roll2c(f, j,1)*a6i[j+4]/(a6i[0]*dx) #j+r->j-r and +->-
        #v=(bondary_roll2c(f,-1,1)-f)/dx
    return v

def U_RK3(U,dx,dt,t,gamma):
    '''
    U'(x) is given before by function deriv_x_p;
    U_next - U = - \partial F(U) / \partial x *dt;
    now calcu in t with 3-order RK
    '''
    npts=t/dt
    #gamma=1.4
    for j in np.arange(npts):
        U0=U
        # 1st stage
        c= c_sound(U, gamma)
        lamda_p, lamda_m= lamda(U,c)
        fp=flow(U,lamda_p,c,gamma)
        fm=flow(U,lamda_m,c,gamma)
        F_x=deriv_x(fp,1,dx)+deriv_x(fm,0,dx)
        U=U0-dt*F_x
        # 2nd stage
        c= c_sound(U, gamma)
        lamda_p, lamda_m= lamda(U,c)
        fp=flow(U,lamda_p,c,gamma)
        fm=flow(U,lamda_m,c,gamma)
        F_x=deriv_x(fp,1,dx)+deriv_x(fm,0,dx)
        U=0.75*U0+0.25*(U-dt*F_x)

        # 3rd stage
        c= c_sound(U, gamma)
        lamda_p, lamda_m= lamda(U,c)
        fp=flow(U,lamda_p,c,gamma)
        fm=flow(U,lamda_m,c,gamma)
        F_x=deriv_x(fp,1,dx)+deriv_x(fm,0,dx)
        U=(U0+2*(U-dt*F_x))/3
    return U

def l2error_each(U1,U2):
    n=len(U1)
    U_l2er=np.zeros((n,1))
    for i in np.arange(n):
        U_l2er[i]=np.sqrt( sum( (U1[i]-U2[i])**2 )/n )
    return U_l2er

################## main execution function 
if __name__ == '__main__':
    ########## x,t mesh:
    domain=(-0.5, 0.5, .0)
    dx = 1./npxs
    x=np.arange(domain[0],domain[1],dx)

    t1=0.14 #t0=0

    ########## initial:
    left0  = (1.,  1.,    0.) #but here the order is [p, rho, u]!!!!
    right0 = (0.1, 0.125, 0.)

    U  =np.zeros([3,npxs]) # U=[rho, u, p]^T at t=0
    left1=np.roll(left0,-1); right1=np.roll(right0,-1);
    for i in np.arange(len(U)):
        U[i,x<=0]=left1[i]
        U[i,x> 0]=right1[i]

##################### solve U at t=t1, the exact and the fvs:
    #### exact, U1e
    positions, regions, values = exact.solve( left_state=left0,  right_state=right0, geometry=domain, t=t1, gamma=gamma, npts=npxs)
    U1e=np.zeros([3,npxs])
    U1e[0]=values['rho']; U1e[1]=values['u']; U1e[2]=values['p'];
    #### fvs, U1
    U1=np.zeros([3,npxs])
    U1=U_RK3(U,dx,dt,t1,gamma)
    print('U1e is \n',U1e[2],'\n U1 is \n', U1[2])

################### plot
    p1 = plt.subplot(2, 2, 1)
    plt.plot(x,U1e[0],label='U1e',linewidth=2)
    plt.plot(x,U1 [0],label='U1',color='r',linewidth=1.5)
    plt.legend()
    plt.xlabel('x') 
    plt.ylabel(['rho(t1)=',t1])
    plt.title('density-x')
    
    p2 = plt.subplot(2, 2, 2)#plt.figure(2)
    plt.plot(x,U1e[1],label='U1e',linewidth=2)
    plt.plot(x,U1 [1],label='U1',color='r',linewidth=1.5)
    plt.legend()
    plt.xlabel('x') 
    plt.ylabel(['u(t1)=',t1])
    plt.title('velocity-x')

    p3 = plt.subplot(2, 2, 3)
    plt.plot(x,U1e[2],label='U1e',linewidth=2)
    plt.plot(x,U1 [2],label='U1',color='r',linewidth=1.5)
    plt.legend()
    plt.xlabel('x') 
    plt.ylabel(['p(t1)=',t1])
    plt.title('pressure-x')
    
    plt.show()

################## print l2-mod error:
    print(l2error_each(U1,U))
# end prog
