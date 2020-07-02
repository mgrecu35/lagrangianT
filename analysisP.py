import glob
from netCDF4 import Dataset
from numpy import *
files=sorted(glob.glob('trajQ*.4'))
import julia
jl=julia.Julia()
jl.include("gridPrecip.jl")
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

fERA=Dataset('evap_precip201511.nc','r')
fIMERG=Dataset('precipNov11_2015.nc','r')
xi=fIMERG['lon'][:]
yi=fIMERG['lat'][:]
ximerg=[]
yimerg=[]
for i in range(600):
    ximerg.append(xi[i*6:i*6+6].mean())
for i in range(300):
    yimerg.append(yi[i*6:i*6+6].mean())
p_imerg=fIMERG['precip'][:,:,:]

p_era=fERA['tp'][:,:,:]
for i in range(1,480,2):
    p_era[i,:,:]=p_era[i,:,:]-p_era[i-1,:,:]

x_era=fERA['longitude'][:]
y_era=fERA['latitude'][:]
from mpl_toolkits.basemap import Basemap, cm
import matplotlib.pyplot as plt

plt.figure()
plt.pcolormesh(p_era[0,::-1,:],norm=matplotlib.colors.LogNorm())
plt.figure()
plt.pcolormesh(p_imerg[0,:,:].T,norm=matplotlib.colors.LogNorm())

pmap_imerg=zeros((360,180))
pmap_era=zeros((360,180))
cmap_imerg=zeros((360,180))
cmap_era=zeros((360,180))
src=zeros((360,90,20,41),float)
csrc=zeros((360,90,20,41),float)
old=zeros((360,90,20),float)
cold=zeros((360,90,20),float)
ic=0
import datetime
for fn in files[:]:
    ic+=1
    xlon,ylat,p,q,f,T,bdate=readf(fn)
    bdateS=bdate.data[0][0,0]
    yy=int(bdateS[0])
    mm=int(bdateS[1])
    dd=int(bdateS[2])
    hh=int(bdateS[3])
    dt=int(bdateS[5])
    a=nonzero(p[0,0,0,:]-p[1,0,0,:]<-30)
    b=nonzero(p[-1,0,0,:][a]>700)
    p1=p[:,0,0,a[0][b]]
    q1=q[:,0,0,a[0][b]]
    x1=xlon[:,0,0,a[0][b]]
    y1=ylat[:,0,0,a[0][b]]
    stime=datetime.datetime(2015,mm,dd,hh)-datetime.timedelta(minutes=-dt)
    dt2=stime-datetime.datetime(2015,11,1,0)
    ih=int(dt2.total_seconds()/3600/6)
    print(ih)
    p_era1d,p_imerg1d,cmap_era,pmap_era,\
    cmap_imerg,pmap_imerg=\
                           jl.gridPrecip(pmap_era,cmap_era,\
                                         pmap_imerg,cmap_imerg,\
                                         x1,y1,p_era,p_imerg,ih)

    #t1=T[:,0,0,a[0][b]]
    #tK1=t1*(p1/1000.)**(0.287)
    #rho1=p1/287./tK1*1e2
    #print(len(b[0]))
    #stop
a=nonzero(array(p_imerg1d)>-1)
print(corrcoef(array(p_imerg1d)[a],array(p_era1d)[a]))
print(array(p_imerg1d)[a].mean()/2.,array(p_era1d)[a].mean()*1000)

m = Basemap(projection='npstere',boundinglat=0,lon_0=0,resolution='l')
from pyresample import geometry, utils, image
import numpy as np
x=arange(361)
y=-90+arange(181)
lons,lats=np.meshgrid(x,y)
x2,y2=m(lons,lats)
a=nonzero(cmap_era>0)
pmap_era[a]/=cmap_era[a]
a=nonzero(cmap_imerg>0)
pmap_imerg[a]/=cmap_imerg[a]

plt.figure(figsize=(8,9))
m.drawcoastlines()
m.drawparallels(np.arange(-80.,81.,20.))
m.drawmeridians(np.arange(-180.,181.,20.))
pmap_eram=np.ma.array(pmap_era,mask=pmap_era<0.1)
pmap_eram=np.ma.array(pmap_eram,mask=cmap_era!=cmap_imerg)
m.pcolormesh(x2,y2,(pmap_eram.T),cmap='jet',vmin=0.1,vmax=5)
plt.title('ERA 12/2- 12/10/2015')
c=plt.colorbar()
c.ax.set_title('mm')
plt.savefig('eraPrecipDec02-10.png')

plt.figure(figsize=(9,8))
m.drawcoastlines()
m.drawparallels(np.arange(-80.,81.,20.))
m.drawmeridians(np.arange(-180.,181.,20.))
pmap_imergm=np.ma.array(pmap_imerg,mask=pmap_imerg<0.1)
pmap_imergm=np.ma.array(0.5*pmap_imergm,mask=cmap_era!=cmap_imerg)
m.pcolormesh(x2,y2,(pmap_imergm.T),cmap='jet',vmin=0.1,vmax=5)
plt.title('IMERG 12/2- 12/10/2015')
c=plt.colorbar()
c.ax.set_title('mm')
plt.savefig('imergPrecipDec02-10.png')

