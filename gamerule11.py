import numpy as np
from vispy import app
from vispy.gloo import clear, set_clear_color, set_viewport, Program
from pylab import *
import matplotlib  
import matplotlib.pyplot as plt

#pre
N=60#numof total agents, ends
l=10 #numof each agt's para
    #assignment for psi(l), by the order of camp[0], life[1], velocity direction1[2],2[3], aera of damage(face)1[4],2[5], position gas1[6],2[7], attack protantial1[10],2[11](no use now, not set), move tendency[8], number[9]. [0][8] constant, [1] initial, [2]~[5] independent, others are dependent
n=20 #space,y,1
m=20 #also,x,2
t=300#time step
life=10#degree of beatproof
psi=zeros(l*n*m)
psi.shape=l,n,m

#main progs#!!rewrite, 6 step in iterate, average field?
def rand123(n):#n belongs integer, return -(n-1)/2:1:n/2 randomly
    a=np.random.random()
    return floor(n*a)

def trsl(A,xd,yd):
    return np.roll(np.roll(A,xd,axis=1),yd,axis=0)

def produce1(N1,N2,l,n,m,life):#produce them randomly
    psi=zeros((l,n,m))
    i=0
    while(i<N1):
        y=rand123(n)
        x=rand123(m)
        if psi[0,y,x]==0:
            #camp=(i-(i/2)*2)*2-1
            psi[:,y,x]=[-1,life, 0,0,0,0, 0,0,0, i]
            i+=1
    i=0
    while(i<N2):
        y=rand123(n)
        x=rand123(m)
        if psi[0,y,x]==0:
            psi[:,y,x]=[1,life, 0,0,0,0, 0,0,0, i]
            i+=1
    return psi

def adjuge(A,B,j,i,ex):
    l,n,m=A.shape
    dj=A[2,j,i]
    di=A[3,j,i]
      #not move if dj*di=0 or someone is here at the last step
    if A[0,j+dj,i+di] or not (dj or di):
        B[:,j,i]=A[:,j,i]
        ex[3]+=1
        return 0
      #not move if another one also want to move here
    a=20/((abs(dj)+abs(di))**2+1)
    #print a,A[:,j+dj,i+di]
    if a>=A[8,j+dj,i+di] and B[0,j+dj,i+di]==0:
        B[:,j+dj,i+di]=A[:,j,i]
        ex[4]+=1
        return 1
    else:
        B[:,j,i]=A[:,j,i]
        ex[5]+=1
        return 0

def iterate(A):#only once
    #assign value
    l,n,m=A.shape
    B=zeros((l,n,m))#only make B empty, because just search Aij which is not empty
    ex=zeros(9)#a good method to adjust prog is set printor everywhere and shrink the area
    A[6:9]=0
      #compu gas_A
    for j in arange(3):
        for i in arange(3):
            if j==1 and i==1:
                continue
            a=trsl(A[0],i-1,j-1)
            a1=a2=a
            a1[a1==-1]=1
            a1[a1!=1]=0
            a2[a2==1]=1
            a2[a2!=1]=0
            A[6]+=a1
            A[7]+=a2
    #decide where for everyone
    for j in arange(n):
        if j+1==n:
            j=-1
        for i in arange(m):
            if i+1==m:
                i=-1
            if A[0,j,i]==0:
                continue
            ex[(A[0,j,i]+1)/2]+=1
            #estimate value and decide, where the changeable subrule is
            if A[(-A[0,j,i]+13)/2,j,i]==0:# or A[0,j,i]==1:#agts in camp=1 goes totally randomly, camp=-1 use new rule!!!!!
                A[2:4,j,i]=[rand123(3)-1,rand123(3)-1]#random move if no enemy around
                continue
            choose=0
            for dj in arange(3):
                if j+dj-1==n:
                    j=1-dj
                for di in arange(3):
                    if i+di-1==m:
                        i=1-di
                    p_we=A[(A[0,j+dj-1,i+di-1]+13)/2,j+dj-1,i+di-1]
                    p_en=A[(-A[0,j+dj-1,i+di-1]+13)/2,j+dj-1,i+di-1]
                    yn1=1-bool(dj*di)
                    yn2=bool(A[0,j+dj-1,i+di-1])
                    yn=yn1*yn2+(1-yn1)*(1-yn2)
                    gs=(100-(p_we-4)**2-(p_en-1)**2)*yn
                    if gs>choose:
                        choose=gs
                        A[2,j,i]=dj-1
                        A[3,j,i]=di-1
    #adjuge moving
      #compu move trendency
    for j in arange(n):
        if j+1==n:
            j=-1
        for i in arange(m):
            if i+1==m:
                i=-1
            a=A[0:4,j,i]
            A[8,j+a[2],i+a[3]]+=10*bool(a[0])/((abs(a[2])+abs(a[3]))**2+1)#diagnal succumb to straight, +1 in denominator to avoid inf, which cannot display in computor-mained but not mathematical signal language, not so good, 10:5:2
      #for eachone, arange it
    for j in arange(n):
        if j+1==n:
            j=-1
        for i in arange(m):
            if i+1==m:
                i=-1
            if A[0,j,i]==0:
                continue
            adjuge(A,B,j,i,ex)
            ex[2]+=1
            ex[6]=ex[2]
    for j in arange(n):
        if j+1==n:
            j=-1
        for i in arange(m):
            if i+1==m:
                i=-1
            if B[0,j,i]==0:
                continue
    #decide face
      #compu gas_B
    for j in arange(3):
        for i in arange(3):
            if j==1 and i==1:
                continue
            a=trsl(B[0],i-1,j-1)
            a1=a2=a
            a1[a1==-1]=1
            a1[a1!=1]=0
            a2[a2==1]=1
            a2[a2!=1]=0
            B[6]+=a1
            B[7]+=a2
      #choose face
    for j in arange(n):
        if j+1==n:
            j=-1
        for i in arange(m):
            if i+1==m:
                i=-1
            if B[0,j,i]==0:
                continue
            for dj in arange(3):
                if j+dj-1==n:
                    j=1-dj
                for di in arange(3):
                    if i+di-1==m:
                        i=1-di
                    if B[0,j+dj-1,i+di-1]==-B[0,j,i]:
                        B[4:6,j,i]=[dj-1,di-1]
                        break#choose the first one, no estimating
    #hurt
    for j in arange(n):
        if j+1==n:
            j=-1
        for i in arange(m):
            if i+1==m:
                i=-1
            if B[0,j,i]==0:
                continue
            for jp in arange(3):
                if j+jp-1==n:
                    j=1-jp
                for ip in arange(3):
                    if i+ip-1==m:
                        i=1-ip
                    a=B[:,j+jp-1,i+ip-1]
                    B[1,j,i]-=bool(a[0])*bool(a[0]-B[0,j,i])*(1-bool(jp-1+a[4]))*(1-bool(ip-1+a[5]))
            if B[1,j,i]<=0:
                ex[(B[0,j,i]+1)/2]-=1
                ex[6]-=1
                B[0:5,j,i]=0
                B[9,j,i]=0
    ex[8]=ex[0]-ex[1]
    #print ex
    if ex[0]*ex[1]<=0:# or ex[2]<N/10:
        ex[7]=1
    return B,ex



# Colormaps
colormaps = np.ones((16, n, 4)).astype(np.float32)
values = np.linspace(0, 1, n)[1:-1]

# Hot colormap
colormaps[0, 0] = 0, 0, 1, 1  # Low values  (< vmin)
colormaps[0, -1] = 0, 1, 0, 1  # High values (> vmax)
colormaps[0, 1:-1, 0] = np.interp(values, [0.00, 0.33, 0.66, 1.00],
                                          [0.00, 1.00, 1.00, 1.00])
colormaps[0, 1:-1, 1] = np.interp(values, [0.00, 0.33, 0.66, 1.00],
                                          [0.00, 0.00, 1.00, 1.00])
colormaps[0, 1:-1, 2] = np.interp(values, [0.00, 0.33, 0.66, 1.00],
                                          [0.00, 0.00, 0.00, 1.00])

# Grey colormap
colormaps[1, 0] = 0, 0, 1, 1  # Low values (< vmin)
colormaps[1, -1] = 0, 1, 0, 1  # High values (> vmax)
colormaps[1, 1:-1, 0] = np.interp(values, [0.00, 1.00],
                                          [0.00, 1.00])
colormaps[1, 1:-1, 1] = np.interp(values, [0.00, 1.00],
                                          [0.00, 1.00])
colormaps[1, 1:-1, 2] = np.interp(values, [0.00, 1.00],
                                          [0.00, 1.00])
# Jet colormap
# ...


img_vertex = """
attribute vec2 position;
attribute vec2 texcoord;

varying vec2 v_texcoord;
void main()
{
    gl_Position = vec4(position, 0.0, 1.0 );
    v_texcoord = texcoord;
}
"""

img_fragment = """
uniform float vmin;
uniform float vmax;
uniform float cmap;

uniform sampler2D image;
uniform sampler2D colormaps;
uniform vec2 colormaps_shape;

varying vec2 v_texcoord;
void main()
{
    float value = texture2D(image, v_texcoord).r;
    float index = (cmap+0.5) / colormaps_shape.y;

    if( value < vmin ) {
        gl_FragColor = texture2D(colormaps, vec2(0.0,index));
    } else if( value > vmax ) {
        gl_FragColor = texture2D(colormaps, vec2(1.0,index));
    } else {
        value = (value-vmin)/(vmax-vmin);
        value = 1.0/512.0 + 510.0/512.0*value;
        gl_FragColor = texture2D(colormaps, vec2(value,index));
    }
}
"""


class Canvas(app.Canvas):
    def __init__(self):
        app.Canvas.__init__(self, size=(512, 512),
                            keys='interactive')
        self.image = Program(img_vertex, img_fragment, 4)
        self.image['position'] = (-1, -1), (-1, +1), (+1, -1), (+1, +1)
        self.image['texcoord'] = (0, 0), (0, +1), (+1, 0), (+1, +1)
        #self.image['vmin'] = +0.1
        #self.image['vmax'] = +0.9
        self.image['vmin'] = 0.1
        self.image['vmax'] = 1.0
        self.image['cmap'] = 1  # Colormap index to use

        self.image['colormaps'] = colormaps
        self.image['colormaps'].interpolation = 'linear'
        self.image['colormaps_shape'] = colormaps.shape[1], colormaps.shape[0]

        self.image['image'] = I.astype('float32')
        self.image['image'].interpolation = 'linear'

        set_clear_color('black')

        self.timer = app.Timer(0.25, self.on_timer)
        self.timer.start()

        self.show()

    def on_resize(self, event):
        width, height = event.physical_size
        set_viewport(0, 0, *event.physical_size)

    def on_draw(self, event):
        clear(color=True, depth=True)
        self.image.draw('triangle_strip')

    def on_timer(self, event):
        #self.clock += 0.001 * 1000.0 / 60.
        #self.program['theta'] = self.clock
        global I
        psi=iterate(psi)
        I=abs(psi[0])
        self.I = I
        #print psi
        self.image['image'] = I.astype('float32')
        self.update()
        #self.timer.stop()

    def on_key_press(self, ev):
        if ev.key.name == 'Space':
            if self.timer.running:
                self.timer.stop()
            else:
                self.timer.start()



#testing
'''
#produce them
n=20
m=20
t=200
N1=30
N2=30
life=10
psi=produce1(N1,N2,l,n,m,life)
#psi[:,3,3]=[2,5,1,-1]#revise directly
#to display
psi1=zeros((t+1,l,n,m))
psi2=zeros((t+1,9))
psi1[0]=psi
psi2[0]=[N1,N2,0,0,0,0,N1+N2,0/1,N1-N2]
tm=t
for i in arange(t):
    itr=iterate(psi)
    psi=itr[0]
    psi1[i+1]=psi
    psi2[i+1]=itr[1]
    if psi2[i+1,7]==1:
        tm=i+1
        break


#display
  #in cmd
print 'camp in every step is: ','\n'
print tm
#print psi1[:,0:2,:,:]
  #in plot
p0=plt.plot(arange(tm),psi2[0:tm,6],'g')#,psi2[:,0],psi2[:,1])
hold
p1=plt.plot(arange(tm),psi2[0:tm,0],'r*')
hold
p2=plt.plot(arange(tm),psi2[0:tm,1],'b*')
hold
p3=plt.plot(arange(tm),psi2[0:tm,8],'g')
plt.xlabel('time step: t')
plt.ylabel('number: (green add and subt, red c-1, blue c1)')
plt.title('')
plt.axis()
y_min=0
if min(psi2[8])<y_min:
    y_min=min(psi2[8])
plt.ylim(y_min-1,N+1)
plt.show()
'''

  #in Canvas
'''psi=psi1[0]
I=psi[0]
if __name__ == '__main__':
    canvas = Canvas()
    app.run()'''

#statis
tn=10
csqc=zeros((tn,2))
for j in arange(tn):
    n=10
    m=10
    t=200
    N1=30
    N2=30
    life=10
    psi=produce1(N1,N2,l,n,m,life)
    psi1=zeros((t+1,l,n,m))
    psi2=zeros((t+1,9))
    psi1[0]=psi
    psi2[0]=[N1,N2,0,0,0,0,N1+N2,0/1,N1-N2]
    for i in arange(t):
        itr=iterate(psi)
        psi=itr[0]
        psi1[i+1]=psi
        psi2[i+1]=itr[1]
        if psi2[i+1,7]==1:
            print psi2[i+1]#####
            csqc[j]=psi2[i+1,0:2]
            break
        csqc[j]=psi2[i+1,0:2]
    print j,csqc[j]
print csqc
cy=(csqc[:,0]-csqc[:,1]+0.0)/(N1+N2+0.0)
print cy
win=0
for i in arange(tn):
    win+=sign(cy[i])
print win
plt.plot(arange(tn),cy,'*')
hold
pzreo=plt.plot([0,tn],[0,0],'r')
plt.xlabel('time step: t')
plt.ylabel('the final win-index of c1 in time: ((N1-N2)/(N1+N2))_tm')
plt.ylim(-1.1,1.1)
plt.xlim(-0.1,tn-0.9)
show()

