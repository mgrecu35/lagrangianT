import glob
from netCDF4 import Dataset
from numpy import *
files=sorted(glob.glob('trajQ*.4'))
import julia
jl=julia.Julia()
jl.include("sourceAttrib.jl")
def readf(fn):
    f=Dataset(fn,'r')
    x=f['lon'][:,:,:]
    y=f['lat'][:,:,:]
    p=f['p'][:,:,:]
    q=f['Q'][:,:,:]
    T=f['TH'][:,:,:]
    bdate=f['BASEDATE'][:,:,:,:]
    return x,y,p,q,f,T,bdate 

import matplotlib
matplotlib.rcParams.update({'font.size': 14})


from mpl_toolkits.basemap import Basemap, cm
import matplotlib.pyplot as plt
cmap=zeros((360,120))
cmap2=zeros((360,120))
src=zeros((360,90,20,41),float)
csrc=zeros((360,90,20,41),float)
old=zeros((360,90,20),float)
cold=zeros((360,90,20),float)
ic=0
oldS=0.
evapS=0.
for fn in files[0:]:
    ic+=1
    xlon,ylat,p,q,f,T,bdate=readf(fn)
    a=nonzero(p[0,0,0,:]-p[1,0,0,:]<-30)
    b=nonzero(p[-1,0,0,:][a]>700)
    p1=p[:,0,0,a[0][b]]
    q1=q[:,0,0,a[0][b]]
    x1=xlon[:,0,0,a[0][b]]
    y1=ylat[:,0,0,a[0][b]]
    t1=T[:,0,0,a[0][b]]
    tK1=t1*(p1/1000.)**(0.287)
    rho1=p1/287./tK1*1e2
    print(len(b[0]))
    dp=50.
    src,csrc,old,cold,oldM,evapM=jl.attrib(p1,rho1,x1,y1,tK1,q1,dp,src,csrc,old,cold)
    oldS+=oldM
    evapS+=evapM
    print(oldS/(oldS+evapS))
    #continue
    for i in a[0][b]:
        ix=int(xlon[0,0,0,i]+180)
        iy=int(ylat[0,0,0,i]+30)
        cmap[ix,iy]+=1
        ix2=int(xlon[-1,0,0,i]+180)
        iy2=int(ylat[-1,0,0,i]+30)
        if ix2>=0 and iy2>=0:
            cmap2[ix2,iy2]+=1
import matplotlib.pyplot as plt
from numpy import *
import numpy as np

def onclick(event):
    if not event.dblclick:
        pos.append([event.xdata,event.ydata])
        #print(pos)
    if event.dblclick:
        pL.append(pos)

        

m = Basemap(projection='npstere',boundinglat=0,lon_0=0,resolution='l')
from pyresample import geometry, utils, image

import matplotlib
iplot=1
for itime in range(41):
    for ilev in range(0,20):
        a=nonzero(csrc[:,:,ilev,itime]>0)
        src[:,:,ilev,itime][a]/=csrc[:,:,ilev,itime][a]


src1=src[:,:,:,:].sum(axis=3)
plt.figure(figsize=(10,10))
m.drawcoastlines()
m.drawparallels(np.arange(-80.,81.,20.))
m.drawmeridians(np.arange(-180.,181.,20.))
x=-179.5+arange(360)
y=0.5+arange(90)
lons,lats=np.meshgrid(x,y)
x2,y2=m(lons,lats)
ic=1
evap=src1[:,:,:].sum(axis=2)*dp/1e3/9.81*1e2/1e3
evapm=ma.array(evap,mask=evap<0.0001)
m.pcolormesh(x2,y2,evapm.T/ic*1e3,cmap='jet',vmax=20.0,\
             norm=matplotlib.colors.LogNorm())
plt.title('Evaporation')
c=plt.colorbar()
c.ax.set_title('mm')
plt.savefig('evaporation.png')
for ilev in range(0,20):
    a=nonzero(cold[:,:,ilev]>0)
    old[:,:,ilev][a]/=cold[:,:,ilev][a]
    
plt.figure(figsize=(10,10))
m.drawcoastlines()
m.drawparallels(np.arange(-80.,81.,20.))
m.drawmeridians(np.arange(-180.,181.,20.))
x=-179.5+arange(360)
y=0.5+arange(90)
lons,lats=np.meshgrid(x,y)
x2,y2=m(lons,lats)


oldM=old[:,:,:].sum(axis=2)*dp/9.81/1e3*1e2/1e3
#oldM=cold[:,:,:11].sum(axis=2)
oldMm=ma.array(oldM,mask=oldM<0.0001)
m.pcolormesh(x2,y2,oldMm.T/ic*1e3,cmap='jet',vmax=20.0,\
             norm=matplotlib.colors.LogNorm())
plt.title('Old moisture')
c=plt.colorbar()
c.ax.set_title('mm')
plt.savefig('oldMoisture.png')
if iplot==1:
    plt.figure(figsize=(11,11))
    m.drawcoastlines()
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))
    x=-179.5+arange(360)
    y=-29.5+arange(120)
    lons,lats=np.meshgrid(x,y)
    x2,y2=m(lons,lats)
    m.pcolormesh(x2,y2,cmap.T,cmap='jet', norm=matplotlib.colors.LogNorm())
    plt.title('Distribution of final locations')
    plt.colorbar()
    plt.savefig('densT0Loc.png')
    plt.figure()
    pos=[]
    pL=[]
    f,a = plt.subplots(1,1,figsize=(11,11))
    f.canvas.mpl_connect('button_press_event', onclick)
    m.drawcoastlines()
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))
    x=-179.5+arange(360)
    #y=0.5+arange(90)
    lons,lats=np.meshgrid(x,y)
    x2,y2=m(lons,lats)
    m.pcolormesh(x2,y2,cmap2.T,cmap='jet', norm=matplotlib.colors.LogNorm())
    plt.title('Distribution of initial locations')
    plt.colorbar()
    plt.savefig('densT0-10daysLoc.png')
