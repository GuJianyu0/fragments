##############python学习与测试
import numpy as np
from numpy import *
import matplotlib.pyplot as plt
#from sympy.abc import ab1,ab2

############# global
a=[[0,1,2,3],[7,7,7,7],[0,-1,-2,-3]]
a=np.array(a)
c=matrix(a)
#print(a,c)
hbar=1

############# functions
def a123(a):
    return hbar*abs(a)

def bondary_roll2c(a,n):
    b=np.roll(a,n,1)
    if n>=0: #translate right n
        for i in np.arange(n):
            b[:,i]=a[:,0]
    if n<0: #translate left |n|
        for i in np.arange(n,0):
            b[:,i]=a[:,-1]
    return b

############# test
a[a<0]=0
print(a)#,'\n',b)


############## tested
#a=[0,1,ab1]
#bool() is just a function in python, not class of var
#print('u1 is \n',u1)
#b=a[-1:-1]
#c=a[-2:-1]
#d=a[-1:-2]
#a[-1:-1]=a[1]
