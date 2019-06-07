#A simple program of plotting C-shape naca0012 square-grid mesh, 
#by solving variable coefficients Poisson equation of coordinates transformation#and by iteration matrix equation from boundary mesh.
#By Jy.
import numpy as np
import matplotlib.pyplot as plt

####instructions and parameters
'''
curve and initial_con, transform/equations, iterations, mesh
网格坐标和真实坐标满足一个非常系数Laplace方程, 然后保持边界网格不变迭代求解这个方程, 求出内部网格

N~j~eita~沿着放射外标度~row, M~i~ksi~沿着长弯C标度~column, 是网格坐标
(X,Y)是真是空间的(x,y)值
离散计算空间, 改变下标: X(eita,ksi)=X[j,i]
'''
M=100 #ksi num i; M,N奇偶影响不大
N=30  #eita num j


####functions
#1.probelms setting...
def naca0012(x):
    return 0.1781*x**0.5 -0.0756*x**2 -0.1705*x**3 +0.0609*x**4 +0.0072  #这个拟合函数有点问题, x=0时翘起来超过x轴了, 把那个差量加回来, 影响不大

def initXY(N,M):
    '''
    用 naca0012机翼 和 C型网格
    '''
    #待求值不要统一用zeros或ones, 否则会出现C_ij分子为0, 用random就好:
    X=(np.random.rand(N,M)-0.5)*16; Y=(np.random.rand(N,M)-0.5)*8
    #X=np.zeros((N,M)); Y=np.zeros((N,M)); #do not use X=Y=..., which is Id X=Y

    #变eita:
    for j in range(0,N):
        Y[j,0]=2.*j/(N-1)  #C型右上边线的真实y坐标
        X[j,0]=4           #C型右上边线的真实x坐标

        Y[j,M-1]=-Y[j,0]   #C型 右下 y
        X[j,M-1]=4         #C型 右上 x

    #变ksi, C型外圈线, 分三段:
    M1=int(M*3/10); M2=int(M-M1); #取整, 这里是M1=15, M2=35; 共M点
    for i in range(0,M1): #外圈上边, i=0~(M1-1)闭, 15点14段
        Y[N-1,i]=2
        X[N-1,i]=4-4.*i/(M1-1 -0) #(i=0, X=4; i=14, X=0;)
    for i in range(M1,M2):#外圈左边, i=M1~(M2-1)闭, 23点22段
        X[N-1,i]=-2*np.sin( np.pi*(i-M1+1)/(M2-M1+1) ) #(i=15, X=sin(pi*1/21); i=34, X=sin(pi*20/21))
        Y[N-1,i]=2*np.cos( np.pi*(i-M1+1)/(M2-M1+1) )
    for i in range(M2,M): #外圈下边, i=M2~(M-1)闭, 15点14段
        Y[N-1,i]=-2
        X[N-1,i]=4.*(i-M2)/(M-1 -M2) #(i=M2, X=0; i=M-1, M=4;)
    
    #变ksi, C型内圈线, 分4段:
    M1=int(M/5); M3=int(M-M/5); M2=int(M/2); #M2=25, 当M=50
    for i in range(0,M1): #机翼右侧区域, 去向往-x, 没够到机翼
        Y[0,i]=0
        X[0,i]=4-3.*i/(M1-1 -0) #*(1-1./M1) #(i=0, X=4; i=M1-1, x=1+small )
        ##print(Y[0,i]," 第1段")
    for i in range(M1,M2): #机翼所在区域, 上 #!!!!
        X[0,i]=1-1.*(i-M1)/(M2-1 -M1) #*(1-1./(M3-M1)) #(i=M1, X=1; i=M2-1, X=1+small)
        Y[0,i]=naca0012(X[0,i]) #naca0012拟合曲线
        ##print(Y[0,i]," 第2段")
    for i in range(M2,M3): #机翼所在区域, 下
        X[0,i]=X[0,M-1 -i] #对称复制, 变换式可由线性方程拟出, 特别地 M=2*M2
        Y[0,i]=-Y[0,M-1 -i]
    for i in range(M3,M): #机翼右侧区域的x坐标, 来向往+x, 没够到机翼
        Y[0,i]=0
        X[0,i]=X[0,M-1 -i]
    
    return X,Y

#2.iteration
#系数 C_ij, D_ij
def moveMat(A,iP,jP):
    '''
    a_{i,j} --> return a_{i+iP,j+jP}
    a very useful funciton
    '''
    A=np.roll(A,-jP,axis=0)
    A=np.roll(A,-iP,axis=1)
    return A

def ABG(X,Y):
    X_ksi =(moveMat(X,1,0)-moveMat(X,-1,0))/2
    Y_ksi =(moveMat(Y,1,0)-moveMat(Y,-1,0))/2
    X_eita=(moveMat(X,0,1)-moveMat(X,0,-1))/2
    Y_eita=(moveMat(Y,0,1)-moveMat(Y,0,-1))/2
    Alpha = X_eita**2+Y_eita**2
    Beta  =-(X_ksi*X_eita+Y_ksi*Y_eita)
    Gamma = X_ksi **2+Y_ksi **2
    return Alpha,Beta,Gamma

#单次迭代式子, 带系数
def JacobiIter1(X,Y):
    Alpha,Beta,Gamma=ABG(X,Y)

    #coefficients C:
    #def coeff(Alpha,Beta,Gamma) #return C
    C=dict={
            (0,0):-2*(Alpha+Gamma),
            (1,0):Alpha, (-1,0):Alpha,
            (0,1):Gamma, (0,-1):Gamma,
            (1,1):Beta/2,(-1,-1):Beta/2,
            (-1,1):-Beta/2,(1,-1):-Beta/2
            }
    
    #Jacobi iteration:
    X_next=-( C[-1,0]*moveMat(X,-1,0) + C[1,0]*moveMat(X,1,0) + C[0,1]*moveMat(X,0,1) + C[0,-1]*moveMat(X,0,-1) )/C[0,0] #小心分母 C_ij==0
    Y_next=-( C[-1,0]*moveMat(Y,-1,0) + C[1,0]*moveMat(Y,1,0) + C[0,1]*moveMat(Y,0,1) + C[0,-1]*moveMat(Y,0,-1) )/C[0,0]
    ##print(X_next)

    #边界不能变(这里整体roll会变了, 就要弄回来)
    X_next[0,:]=X[0,:]; X_next[-1,:]=X[-1,:];
    X_next[:,0]=X[:,0]; X_next[:,-1]=X[:,-1];
    Y_next[0,:]=Y[0,:]; Y_next[-1,:]=Y[-1,:];
    Y_next[:,0]=Y[:,0]; Y_next[:,-1]=Y[:,-1];

    return X_next,Y_next

#K次GS迭代
def IterK(X0,Y0,K):
    X_k,Y_k=X0,Y0 #initial value used is this??
    X_knext,Y_knext=X0,Y0
    for k in range(0,K):
        X_knext,Y_knext=JacobiIter1(X_k,Y_k)
        X_k,Y_k=X_knext,Y_knext
    return X_knext,Y_knext

#3.画出每条线
#一条条地画出网格, 按(ksi,eita)离散编号i,j, 画出其在(X,Y)的线条
def plotmesh(X,Y):
    (N,M)=X.shape
    if ((N,M)!=Y.shape):
        exit()
        return 0
    plt.figure()
    for j in range(0,N):
        plt.plot(X[j,:],Y[j,:])#,'b')
    for i in range(0,M):
        plt.plot(X[:,i],Y[:,i])#,'r')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('naca0012 squaregrid mesh by Laplace solving and iterations')
    plt.show()
    return 0


#main
if __name__ == '__main__':
    #赋初值
    X,Y=initXY(N,M)

    #迭代
    X,Y=IterK(X,Y,2000) #这里Jacobi迭代次数不足10^3效果不好, 看来单步迭代不优秀

    #画图
    plotmesh(X,Y)

#endmain
