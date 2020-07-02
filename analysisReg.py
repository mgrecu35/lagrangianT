from mpl_toolkits.basemap import Basemap, cm
import matplotlib.pyplot as plt
        

m = Basemap(projection='npstere',boundinglat=0,lon_0=0,resolution='l')
from pyresample import geometry, utils, image

from numpy import *
import numpy as np
x=-179.5+arange(360)
y=0.5+arange(90)
lons,lats=np.meshgrid(x,y)
x2,y2=m(lons,lats)
cmap=zeros((90,360),float)
for i in range(90):
    for j in range(360):
        m1=m.is_land(x2[i,j],y2[i,j])
        if m1 is True:
            cmap[i,j]=1

import pickle
pickle.dump(cmap,open('cmap.pklz','wb'))
m.drawcoastlines()
m.drawparallels(np.arange(-80.,81.,20.))
m.drawmeridians(np.arange(-180.,181.,20.))
x=-179.5+arange(360)
y=0.5+arange(90)
lons,lats=np.meshgrid(x,y)
x2,y2=m(lons,lats)
import matplotlib
m.pcolormesh(x2,y2,cmap,cmap='jet', norm=matplotlib.colors.LogNorm())
