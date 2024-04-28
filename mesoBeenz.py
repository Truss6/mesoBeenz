'''

 .----------------.  .----------------.  .----------------.  .----------------.  .----------------.  .----------------.  .----------------.  .-----------------. .----------------. 
| .--------------. || .--------------. || .--------------. || .--------------. || .--------------. || .--------------. || .--------------. || .--------------. || .--------------. |
| | ____    ____ | || |  _________   | || |    _______   | || |     ____     | || |   ______     | || |  _________   | || |  _________   | || | ____  _____  | || |   ________   | |
| ||_   \  /   _|| || | |_   ___  |  | || |   /  ___  |  | || |   .'    `.   | || |  |_   _ \    | || | |_   ___  |  | || | |_   ___  |  | || ||_   \|_   _| | || |  |  __   _|  | |
| |  |   \/   |  | || |   | |_  \_|  | || |  |  (__ \_|  | || |  /  .--.  \  | || |    | |_) |   | || |   | |_  \_|  | || |   | |_  \_|  | || |  |   \ | |   | || |  |_/  / /    | |
| |  | |\  /| |  | || |   |  _|  _   | || |   '.___`-.   | || |  | |    | |  | || |    |  __'.   | || |   |  _|  _   | || |   |  _|  _   | || |  | |\ \| |   | || |     .'.' _   | |
| | _| |_\/_| |_ | || |  _| |___/ |  | || |  |`\____) |  | || |  \  `--'  /  | || |   _| |__) |  | || |  _| |___/ |  | || |  _| |___/ |  | || | _| |_\   |_  | || |   _/ /__/ |  | |
| ||_____||_____|| || | |_________|  | || |  |_______.'  | || |   `.____.'   | || |  |_______/   | || | |_________|  | || | |_________|  | || ||_____|\____| | || |  |________|  | |
| |              | || |              | || |              | || |              | || |              | || |              | || |              | || |              | || |              | |
| '--------------' || '--------------' || '--------------' || '--------------' || '--------------' || '--------------' || '--------------' || '--------------' || '--------------' |
 '----------------'  '----------------'  '----------------'  '----------------'  '----------------'  '----------------'  '----------------'  '----------------'  '----------------' 

'''
##### Imports
import numpy as np
from math import *
import pandas as pd
import os,csv,time,random
import matplotlib.patches as Patches

##### Plotting presets
import matplotlib as mpl
mpl.rcParams['figure.figsize']=[6,4] # set fig size (in inches)
mpl.rcParams['figure.dpi']=300 # set resolution (debugging: 100, draft: 300, final: 900)
mpl.rcParams['font.size']=12 # set global font size
mpl.rcParams['font.family']='times new roman' # set global font
# in a linux OS, this may throw an error, consider changing to 'FreeSerif'
mpl.rcParams['legend.fontsize']=12 # decrease legend font size
mpl.rcParams['lines.markersize']=10 # increase default marker size
mpl.rcParams['mathtext.fontset'] = 'cm' # set math font to computer-modern
from matplotlib import pyplot as plt

trigs=lambda th: [cos(th),sin(th)]
rand_th=lambda:random.random()*2*pi

def cross_2d(a,b):
    return a[0]*b[1]-a[1]*b[0]
def dot_2d(a,b):
    return a[0]*b[0]+a[1]*b[1]

def A_been(r,L):
    if L>2*r: return 2*pi*r**2
    phi=acos(L/(2*r))
    # A_tri=L*R*sin(phi)
    # A_slice=(pi*r**2)*phi/(2*pi)
    # A_crust=2*(A_slice-A_tri)
    # A_circs=2*pi*r**2
    # A_been=A_circs-A_crust
    # # return r*(2*pi*r+L*sin(phi)-phi*r)
    # return A_been
    return 2*pi*r**2 - 2*(phi*r**2-0.5*L*r*sin(phi))
    #2*pi*r**2 -(phi*r**2-L*r*sin(phi))
Acirc=lambda r: pi*r**2

class Domain:
    def __init__(self,x,y):
        self.x=x
        self.y=y
        self.area=(x[1]-x[0])*(y[1]-y[0])
    def rand(obj,N):
        x=np.random.rand(N)*(obj.x[1]-obj.x[0])+obj.x[0]
        y=np.random.rand(N)*(obj.y[1]-obj.y[0])+obj.y[0]
        return x,y
    def update(self,h=1):
        exit_flag=True
        err=0
        Nerrs=0
        N=self.system.N_beenz
        for i in range(N):
            # F=[0,0]
            # M=0
            # x0=self.system.beenz[i].x
            # y0=self.system.beenz[i].y
            # reading coordinates here was less stable
            L0=self.system.beenz[i].L
            # th0=self.system.beenz[i].th
            r0=self.system.beenz[i].r
            A0=2*r0**2#A_been(r0,L0)
            I0=(r0*L0)**2
            for j in range(N):
                x1=self.system.beenz[j].x 
                y1=self.system.beenz[j].y
                L1=self.system.beenz[j].L
                th1=self.system.beenz[j].th
                r1=self.system.beenz[j].r
                # A1=2*r1^2#A_been(r1,L1)
                # I1=(r1*L1)**2
                if i!=j:
                    for si in [+1,-1]:
                        x0=self.system.beenz[i].x
                        y0=self.system.beenz[i].y
                        # L0=self.system.beenz[i].L
                        th0=self.system.beenz[i].th
                        # r0=self.system.beenz[i].r
                        # A0=A_been(r0,L0)
                        # I0=(r0*L0)**2
                        
                        xA=x0+si*L0*cos(th0)*0.5
                        yA=y0+si*L0*sin(th0)*0.5
                        L_0=[t*si*L0/2 for t in trigs(th0)]
                        for sj in [+1,-1]:
                            # x1=self.system.beenz[j].x
                            # y1=self.system.beenz[j].y
                            # L1=self.system.beenz[j].L
                            # th1=self.system.beenz[j].th
                            # r1=self.system.beenz[j].r
                            # A1=A_been(r1,L1)
                            # print(sj)
                            xB=x1+sj*L1*cos(th1)*0.5
                            yB=y1+sj*L1*sin(th1)*0.5
                            
                            D=[xB-xA,yB-yA]
                            D_mag=(D[0]**2+D[1]**2)**0.5
                            d_mag=D_mag-r1-r0
                            # print(f'd={D_mag}')
                            if d_mag<self.system.tol:
                                exit_flag=False
                                err+=-d_mag
                                Nerrs+=1
                                # print(f'{i} {si}, {j} {sj}')
                                # print(f'HERE d={d_mag}')
                                e_d=[Di/D_mag for Di in D] # center to center unit vector
                                # note that vector D is parallel to vector d
                                d=[e_d_i*d_mag for e_d_i in e_d] # vector of max overlap
                                # print(d)
                                
                                L_1=[t*sj*L1/2 for t in trigs(th1)]
                                
                                # lin_0=2*dot_2d(L_0,d)/L0
                                # lin_1=2*dot_2d(L_1,d)/L1
                                A1=2*r1**2;I1=(r1*L1)**2
                                
                                ang_0=2*cross_2d(L_0,d)*I1/(I0+I1)
                                ang_1=2*cross_2d(L_1,d)*I0/(I0+I1)
                                
                                # print(f'{i}: {lin_0},{ang_0}')
                                # print(f'{j}: {lin_1},{ang_1}')
                                
                                self.system.beenz[i].x+=h*d[0]*A1/(A0+A1)
                                self.system.beenz[i].y+=h*d[1]*A1/(A0+A1)
                                self.system.beenz[i].th+=h*ang_0
                                # F[0]+=h*d[0]*A1/(A0+A1)
                                # F[1]+=h*d[1]*A1/(A0+A1)
                                # M+=h*ang_0
                                
                                self.system.beenz[j].x+=h*d[0]*A0/(A0+A1)
                                self.system.beenz[j].y+=-h*d[1]*A0/(A0+A1)
                                self.system.beenz[j].th+=-h*1*ang_1
                        # self.system.beenz[i].x+=F[0]
                        # self.system.beenz[i].y+=F[1]
                        # self.system.beenz[i].th+=M
                        h_BC=h
                        if xA-r0<self.x[0]:
                            d=self.x[0]-(xA-r0)
                            if abs(d)>self.system.tol:
                                exit_flag=False
                                err+=d;Nerrs+=1
                            e_d=[1,0]
                            self.system.beenz[i].x+=h_BC*d#*dot_2d(L_0,e_d)/L0
                            # print('xb')
                        elif xA+r0>self.x[1]:
                            d=self.x[1]-(xA+r0)
                            e_d=[-1,0]
                            if abs(d)>self.system.tol:
                                exit_flag=False
                                err+=d;Nerrs+=1
                                self.system.beenz[i].x+=h_BC*d#*dot_2d(L_0,e_d)/L0
                            # print('xt')
                        if yA-r0<self.y[0]:
                            d=self.y[0]-(yA-r0)
                            e_d=[0,1]
                            if abs(d)>self.system.tol:
                                exit_flag=False
                                err+=d;Nerrs+=1
                                self.system.beenz[i].y+=h_BC*d#*dot_2d(L_0,e_d)/L0                         
                            # print('yb')
                        elif yA+r0>self.y[1]:
                            d=self.y[1]-(yA+r0)
                            e_d=[0,-1]
                            if abs(d)>self.system.tol:
                                exit_flag=False
                                err+=d;Nerrs+=1
                            self.system.beenz[i].y+=h_BC*d#*dot_2d(L_0,e_d)/L0                       
                            # print('yt')
        return exit_flag,err,Nerrs
    # def assign(self)
    def plot(self,i):
        plt.figure()
        ax = plt.gca()
        plt.xlim(self.x)
        plt.ylim(self.y)
        self.system.plot(ax)
        # plt.savefig('beenz'+str(i).zfill(6)+'.png',bbox_inches='tight')
        # plt.close()
class system:
    N_beenz=None
    beenz=[]
    tol=0
    def __init__(self,beenz):
        self.beenz=beenz
        self.N_beenz=len(beenz)
        self.tol=1e-3*min([been.r for been in self.beenz])
    def plot(self,ax):
        for been in self.beenz:
            been.plot(ax)
    def write(self,fname):
        print("Wrtiting "+fname)
        with open(fname,'w') as out:
            out.write('x,y,r,id\n')
            for been in self.beenz:
                for s in [+1,-1]:
                    line=str(been.x+s*0.5*been.L*cos(been.th))+','+str(been.y+s*0.5*been.L*sin(been.th))+','+str(been.r)+','+str(been.i)+'\n'
                    out.write(line)
class been:
    x=None
    y=None
    r=None
    L=None
    i=None
    c=None
    def __init__(self,x,y,r,L,th,i,c=None):
        self.x=x
        self.y=y
        self.r=r
        self.i=i
        self.L=L
        self.th=th
        if c==None:
            self.c=np.random.rand(3)
        else:
            self.c=c
    def plot(self,ax):
        for s in [-1,+1]:
            if s==1: C=self.c
            else: C=[0.5+clr/2 for clr in self.c]
            circle=plt.Circle((self.x+s*0.5*self.L*cos(self.th),self.y+s*0.5*self.L*sin(self.th)),self.r,color=C,alpha=1,ec=self.c,fc=C,lw=2,ls='-')
            ax.add_patch(circle)
def makeBeenz(x,y,r,L,th,i=None,c=None): # this initializes a set of points
    beenz=[]
    if len(x)!=len(y):
        print('ERROR1')
        return None
    else:
        N_beenz=len(x)
    if len(r)==1:
        r=[r for _ in range(N_beenz)]
    if c==None:
        c=[None for _ in range(N_beenz)]
    for X,Y,R,l,TH,I,C in zip(x,y,r,L,th,i,c):
        beenz.append(been(X,Y,R,l,TH,I,C)) 
    return beenz
'''
 .----------------.  .----------------.  .----------------.  .-----------------.
| .--------------. || .--------------. || .--------------. || .--------------. |
| | ____    ____ | || |      __      | || |     _____    | || | ____  _____  | |
| ||_   \  /   _|| || |     /  \     | || |    |_   _|   | || ||_   \|_   _| | |
| |  |   \/   |  | || |    / /\ \    | || |      | |     | || |  |   \ | |   | |
| |  | |\  /| |  | || |   / ____ \   | || |      | |     | || |  | |\ \| |   | |
| | _| |_\/_| |_ | || | _/ /    \ \_ | || |     _| |_    | || | _| |_\   |_  | |
| ||_____||_____|| || ||____|  |____|| || |    |_____|   | || ||_____|\____| | |
| |              | || |              | || |              | || |              | |
| '--------------' || '--------------' || '--------------' || '--------------' |
 '----------------'  '----------------'  '----------------'  '----------------' 
'''
if __name__=='__main__':
    vf=0.75
    D=Domain([-5,5],[-5,5])
    A=D.area*vf
    # x=[]
    # y=[]
    r=[]
    L=[]    
    while A>0:
        R=random.random()*0.5
        l=R*random.random()*2
        r.append(R)
        L.append(l)
        A-=A_been(R,l)
    N=len(r)
    # N=100
    x,y=D.rand(N)
    # x=[-4,4,0]
    # y=[-0.5*cos(pi/4),2,-4]
    
    # r=[random.random()*0.5+0.0 for _ in range(N)]
    # r[1]=1.5
    # L=[1*R+0*random.random()*R*2 for R in r]
    th=[rand_th() for _ in range(N)]
    ids=[i for i in range(N)]
    clr=[]
    for i in range(N):
        c=i/(N-1)
        C=[c**2,4*c*(1-c),(1-c)**2]
        clr.append(C)
    beenz=makeBeenz(x,y,r,L,th,ids,clr)
    D.system=system(beenz)
    D.plot(0)
    err0=inf;h=1;err_min=inf;
    f=plt.figure()
    for i in range(300):
        exit_flag,err,Nerr=D.update(h=1)
        print(f'ERROR({i}): # {Nerr} = {err}')
        # D.plot(i+1)
        if exit_flag:
            print(f'Converged After {i+1} Iterations')
            break
        if err>err0:
            # if solution is diverging, reduce step size
            h=h*(err0/err)**0.5
            print(f'H reduced to: {h}')
        if err<err_min: 
            # only write the outputfile if
            # the current solution is better than the best thus far
            D.system.write('meso.csv')
            err_min=err
        err0=err
        f.gca().plot(i,err,'ko')
    print(f'Remaining Error: {err}')
    plt.xlabel('Cycle (N)')
    plt.ylabel('Overlap (err) [distance]')
    f.savefig('convergence.png',bbox_inches='tight')
    D.plot(i+1)
    
    # plt.xlim([-2,5])
    # plt.ylim([-2,2])