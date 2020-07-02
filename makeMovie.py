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

from mpl_toolkits.basemap import Basemap, cm
import matplotlib.pyplot as plt
cmap=zeros((360,90))
cmap2=zeros((360,90))
src=zeros((360,90,20),float)
csrc=zeros((360,90,20),float)
old=zeros((360,90,20),float)
cold=zeros((360,90,20),float)
for fn in files[9:10]:
    xlon,ylat,p,q,f,T,bdate=readf(fn)
    a=nonzero(p[0,0,0,:]-p[1,0,0,:]<-30)
    b=nonzero(p[-1,0,0,:][a]>700)
    p1=p[:,0,0,a[0][b]]
    q1=q[:,0,0,a[0][b]]
    x1=xlon[:,0,0,a[0][b]]
    y1=ylat[:,0,0,a[0][b]]
    t1=T[:,0,0,a[0][b]]
   
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
iplot=0
it=0
for ilev in range(39,4,-1):
    plt.figure(figsize=(10,10))
    m.drawcoastlines()
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))
    x=-179.5+arange(360)
    y=0.5+arange(90)
    lons,lats=np.meshgrid(x,y)
    x2,y2=m(lons,lats)
    czero=zeros((360,90))
    for x11,y11 in zip(x1[ilev,:],y1[ilev,:]):
        ix=int(x11+180)
        iy=int(y11)
        if ix>=0 and ix<360 and iy>=0 and iy<90:
            czero[ix,iy]+=1

    m.pcolormesh(x2,y2,czero.T,cmap='jet',
                 norm=matplotlib.colors.LogNorm())
    plt.savefig('map%2.2i.png'%it)
    it+=1
