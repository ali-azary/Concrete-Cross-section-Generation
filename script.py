from numpy import *
from matplotlib.pyplot import *
from random import *
import math
import matplotlib.patches as patches
import sys
from subprocess import call
import time

start_time = time.time()

h=1. # mesh size
p=40. # aggregate percentage

# dimensions
lx=50.
ly=50.

# number of vertices of particles
nv=10

# required aggregate area
req=lx*ly*p/100.

f=open('aggregates','w')

# aggregates
fig = figure()
ax1 = fig.add_subplot(111, aspect='equal')
ax1.add_patch(
    patches.Rectangle(
        (0., 0.),   # (x,y)
        lx,          # width
        ly,          # height
        fill=False 
    )
)
centers=[]
rs=[]

theta=linspace(0.,2.*pi,nv+1)
ni=0

# grading
sieves=[1.18,2.36,4.75,9.5,12.5,14.]#19.]
minp=[0.,0.,40.,90.,100.]
maxp=[5.,15.,70.,100.,100.]

# well graded
# cumulative
wellc=[]
for i in range(size(sieves)-1):
    value=minp[i]+(maxp[i]-minp[i])*random()
    if i==0: wellc.append(value)
    if i>0:        
        while value<=wellc[i-1]:
            value=minp[i]+(maxp[i]-minp[i])*random()
        wellc.append(value)
# percentage    
wellp=[wellc[0]]
for i in range(size(sieves)-2):
    wellp.append(wellc[i+1]-wellc[i])

sieves=sieves[::-1]
wellp=wellp[::-1]

for i in range(size(wellp)):
    Ar=wellp[i]/100.*req
    minr=sieves[i+1]/2.
    maxr=sieves[i]/2.
    # generating and placing aggregates
    A=0.
    while sum(A)<Ar:
        elapsed_time = time.time() - start_time
        if elapsed_time > 10:  # If the loop takes more than a minute to run
            call(["python", __file__])  # Restart the script
            exit()  # Exit the current instance of the script
        xc=lx*random()
        yc=ly*random()
        A0=minr+(maxr-minr)*random()
        m=int(4.+6.*random())
        Aj=[]
        alphaj=[]
        p=.9
        b=1.9
        for j in range(m):
            j=j+1
            Aj.append(A0*exp(-p*log(j)-b))
            alphaj.append(2.*pi*random())
        r=A0    
        for j in range(m):
            r=r+Aj[j]*cos(j*theta+alphaj[j])        
        add=1
        x=[]
        y=[]
        tmp=0.
        for i in range(nv):
            tmp=tmp+r[i]*r[i+1]*sin(2.*pi/nv)/2.
        for i in range(nv+1):
            xi=xc+r[i]*cos(theta[i]);
            if xi<=1. or xi>=lx-1.: add=0; break
            yi=yc+r[i]*sin(theta[i]);
            if yi<=1. or yi>=ly-1.: add=0; break
            x.append(xi)
            y.append(yi)    
        if add==1:
            if ni>0:
                for j in range(shape(centers)[0]):
                    ci=centers[j][:]
                    ri=rs[j][:]            
                    di=sqrt((xc-ci[0])**2+(yc-ci[1])**2)
                    if di<max(ri)+max(r): add=0; break
                   
                if add==1:
                    centers.append([xc,yc]); rs.append(r); plot(x,y,color='black'); ni=ni+1
                    A=A+tmp
            if ni==0:
                centers.append([xc,yc]); rs.append(r); plot(x,y,color='black'); ni=ni+1
                A=A+tmp

        if add==1:
            for i in range(nv):
                f.write(str(x[i])+' ')
            for i in range(nv):
                f.write(str(y[i])+' ')
            f.write('\n')
            sys.stdout.write('%i particles with %f percent required area\n'%(ni,A/Ar*100.))
            sys.stdout.flush()

    axis('equal')
    axis('off')
    fig.savefig('section.jpg', dpi=100, bbox_inches='tight')
    
    

f.close()

# writing input file for gmsh for mesh generation
g=open('aggregates','r')
f=open('model.geo','w')

edge1=[]
edge2=[]

f.write('Point(1)={0.,0.,0.,%f};\n'%(h))
f.write('Point(2)={%f,0.,0.,%f};\n'%(lx,h))
f.write('Point(3)={%f,%f,0.,%f};\n'%(lx,ly,h))
f.write('Point(4)={0.,%f,0.,%f};\n'%(ly,h))

f.write('Line(5)={1,2};\n')
f.write('Line(6)={2,3};\n')
f.write('Line(7)={3,4};\n')
f.write('Line(8)={4,1};\n')

edge1.extend([5,6,7,8])

cnt=9
for line in g:
    x=[float(i) for i in real(line.split())[:nv]]
    y=[float(i) for i in real(line.split())[nv:2*nv]]
    for i in range(nv):
        f.write('Point(%i)={%f,%f,0.,%f};\n'%(cnt+i,x[i],y[i],h))    
    for i in range(nv-1):
        f.write('Line(%i)={%i,%i};\n'%(cnt+nv+i,cnt+i,cnt+i+1))
        edge1.append(-(cnt+nv+i))        
    f.write('Line(%i)={%i,%i};\n'%(cnt+2*nv-1,cnt+nv-1,cnt))
    edge1.append(-(cnt+2*nv-1))
    edge2.append([i for i in arange(cnt+nv,cnt+2*nv,1)])
    cnt=cnt+2*nv

f.write('Line Loop(%i)={'%cnt)
for i in range(size(edge1)-1):
    f.write('%i,'%edge1[i])
f.write('%i};\n'%edge1[-1])
f.write('Plane Surface(%i)={%i};\n'%(cnt+1,cnt))
cnt=cnt+2

f.close()    

call(['gmsh','-2','model.geo','-o','matrix.msh'])

f=open('model.geo','a+')

for e in edge2:
    f.write('Line Loop(%i)={'%cnt)
    for i in range(size(e)-1):
        f.write('%i,'%e[i])
    f.write('%i};\n'%e[-1])
    f.write('Plane Surface(%i)={%i};\n'%(cnt+1,cnt))
    cnt=cnt+2

f.close()

# mesh generation using gmsh
call(['gmsh','-2','model.geo','-o','section.msh'])
