#!/usr/bin/env python
# coding: utf-8

# In[64]:


import numpy as np 
from matplotlib import pyplot as plt
import fieldLineTopology as flt
from streamtracer import StreamTracer, VectorGrid
from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import simps
#import waveletRoutines as wr

import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.io import netcdf_file


# In[65]:


# read in file
path = '/home/grads/trcn27/Documents/postdoc/lare3d_jet/Data/'#'C://Chris Research/activeRegion/'
file2read = netcdf_file(path+'0000.nc','r')


# In[66]:


file2read.variables


# In[67]:


bxOg =file2read.variables['bx'][:]*1
byOg =file2read.variables['by'][:]*1
bzOg =file2read.variables['bz'][:]*1


# In[69]:


# smush to centred grid
bx =0.5*(bxOg[:,:,:(bxOg.shape[2]-1)]+bxOg[:,:,1::bxOg.shape[2]])
by =0.5*(byOg[:,:(byOg.shape[1]-1)]+byOg[:,1::byOg.shape[1],:])
bz =0.5*(bzOg[:(bzOg.shape[0]-1),:,:]+bzOg[1::bzOg.shape[0],:,:])

ncells = int(np.size(bx)**(1/3)) + 1
print(ncells)
#Reshape so its 0= x comp 1 =y comp 2 = z comp 
bx = bx.reshape([ncells,ncells,ncells]).transpose(2,1,0)
by = by.reshape([ncells,ncells,ncells]).transpose(2,1,0)
bz = bz.reshape([ncells,ncells,ncells]).transpose(2,1,0)
xv = np.linspace(-10,10,ncells)
yv = np.linspace(-10,10,ncells)
zv = np.linspace(0,20,ncells)
X, Y = np.meshgrid(xv, yv, indexing='ij')
# Flatten the grid arrays to form the input to the interpolator
points = np.vstack([X.ravel(), Y.ravel()]).T
dx = xv[1]-xv[0]
dy = yv[1]-yv[0]
dz = zv[1]-zv[0]
dA = dx*dy
grid_spacing = [dx,dy,dz]


# In[70]:


bField = flt.createSingleField(bx,by,bz)
#bFieldTest = flt.createSingleField(bxRot,byRot,bzRot)


# In[71]:


#get curl
curlField= flt.curl(bField,grid_spacing)

#calculate the winding gauge 
AField = flt.getAFastSingle(bField)

#calculate the winding gauge for the unit speed field 
usf = flt.unitSpeedField(bField,0.01)
AUnit= flt.addDivergenceCleaningTerm(usf,grid_spacing)

AWind = flt.getAFastSingle(AUnit)

# add constant component of A as there be net flux

bzConst = np.sum(bz[:,:,0])/(dA*(ncells-1)*(ncells-1))
AConst = flt.AConst(bzConst,points,dA,[ncells, ncells, ncells])
AField = AField + AConst 



# In[72]:


# calculate the field line helcity and winding densities

testFLHDen =flt.getFLHDenSingle(bField,AField)
testWindDen =flt.getFLHDenSingle(bField,AField)
# twist density
twistDensity = flt.twistDen(bField,curlField,0.1)


# In[73]:


#Check the field comps

fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 4))
axes[0].imshow(bField[:,:,0,0], origin='lower')
axes[1].imshow(bField[:,:,0,1], origin='lower')
axes[2].imshow(bField[:,:,0,2], origin='lower')
fig.tight_layout()
plt.show()


# In[1]:


#Check the curl comps

fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 4))
axes[0].imshow(twistDensity[:,:,0], origin='lower')
axes[1].imshow(twistDensity[:,:,20], origin='lower')
axes[2].imshow(twistDensity[:,:,40], origin='lower')
fig.tight_layout()
plt.show()


# In[61]:


# FLH density, pretty !

fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 4))
axes[0].imshow(testFLHDen[:,:,0].T, origin='lower')
axes[1].imshow(testFLHDen[:,:,20].T, origin='lower')
axes[2].imshow(testFLHDen[:,:,40].T, origin='lower')
fig.tight_layout()
plt.show()


# In[62]:


# FLwind density, pretty !

fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 4))
axes[0].imshow(testWindDen[:,:,10].T, origin='lower')
axes[1].imshow(testWindDen[:,:,20].T, origin='lower')
axes[2].imshow(testWindDen[:,:,40].T, origin='lower')
fig.tight_layout()
plt.show()


# Set up the field line tracer

# In[63]:


grid_spacing = [20/ncells,20/ncells,20/ncells]

#set domain lengths

lx = grid_spacing[0]*ncells
ly = grid_spacing[1]*ncells
lz = grid_spacing[2]*ncells

#set number of points used to calculate the distribution from
nx = 200
ny = 200

# set a minimum strength of field line cut off 
bCut = 0.01

# set z value for anchoring plane (use photosphere z index = 116 here)
zv = 0


# In[ ]:


BxInterp,ByInterp,BzInterp = flt.getInterpolatedFieldSingle(bField,grid_spacing[0],grid_spacing[1],grid_spacing[2])


# In[ ]:


fieldLinesList,goodSeeds,seeds = flt.prepareCurves(bField,BxInterp,ByInterp,BzInterp,grid_spacing,[6,14],[6,14],zv,nx,ny,bCut)


# In[ ]:


flhInterp = flt.getInterpolatedQuantity(testFLHDen,grid_spacing)
flhBzInterp = flt.getInterpolatedQuantity(testFLHDen*bz[:,:,0],grid_spacing)
flwindInterp = flt.getInterpolatedQuantity(testWindDen,grid_spacing)
twistInterp = flt.getInterpolatedQuantity(twistDensity,grid_spacing)


# In[ ]:


flh,indexesflh = flt.fieldLineIntegratedQuantity(flhInterp,goodSeeds,fieldLinesList,seeds,nx,ny)
flhBz,indexesflh = flt.fieldLineIntegratedQuantity(flhBzInterp,goodSeeds,fieldLinesList,seeds,nx,ny)
flw,indexesflw = flt.fieldLineIntegratedQuantity(flwindInterp,goodSeeds,fieldLinesList,seeds,nx,ny)
twistF,indexestwist = flt.fieldLineIntegratedQuantity(twistInterp,goodSeeds,fieldLinesList,seeds,nx,ny)


# In[44]:


plt.imshow(flh,origin='lower',cmap='seismic')
plt.show()


# In[45]:


plt.imshow(flw,origin='lower',cmap='seismic')
plt.show()


# In[46]:


plt.imshow(flhBz,origin='lower',cmap='seismic')
plt.show()


# In[ ]:


plt.imshow(twistF,origin='lower',cmap='seismic')
plt.show()


# In[ ]:


def getFieldLineMetrics(fileLoc,dimensions,grid_spacing):
    #open the file
    file2read = netcdf.NetCDFFile(fileLoc,'r')
    # read the field variables
    bxOg =file2read.variables['bx'][:]*1
    byOg =file2read.variables['by'][:]*1
    bzOg =file2read.variables['bz'][:]*1
    # smush to centred grid
    bx =0.5*(bxOg[:,:,:(bxOg.shape[2]-1)]+bxOg[:,:,1::bxOg.shape[2]])
    by =0.5*(byOg[:,:(byOg.shape[1]-1)]+byOg[:,1::byOg.shape[1],:])
    bz =0.5*(bzOg[:(bzOg.shape[0]-1),:,:]+bzOg[1::bzOg.shape[0],:,:])

    #Reshape so its 0= x comp 1 =y comp 2 = z comp 
    bx = bx.reshape(dimensions).transpose(2,1,0)
    by = by.reshape(dimensions).transpose(2,1,0)
    bz = bz.reshape(dimensions).transpose(2,1,0)

    #spatial grids    
    xv = np.linspace(-10,10,192)
    yv = np.linspace(-10,10,192)
    zv = np.linspace(0,20,192)
    X, Y = np.meshgrid(xv, yv, indexing='ij')
    # Flatten the grid arrays to form the input to the interpolator
    points = np.vstack([X.ravel(), Y.ravel()]).T
    dx = xv[1]-xv[0]
    dy = yv[1]-yv[0]
    dz = zv[1]-zv[0]
    dA = dx*dy
    grid_spacing = [dx,dy,dz]
    
    bField = flt.createSingleField(bx,by,bz)
    
    #get curl
    curlField= flt.curl(bField,grid_spacing)

    #calculate the winding gauge 
    AField = flt.getAFastSingle(bField)

    #calculate the winding gauge for the unit speed field 
    usf = flt.unitSpeedField(bField,0.01)
    AUnit= flt.addDivergenceCleaningTerm(usf,grid_spacing)

    AWind = flt.getAFastSingle(AUnit)

    #add constant component of A as there be net flux

    bzConst = np.sum(bz[:,:,0])/(dA*191*191)
    AConst = flt.AConst(bzConst,points,dA,[192,192,192])
    AField = AField + AConst 
    
    # calculate the field line helcity and winding densities

    testFLHDen =flt.getFLHDenSingle(bField,AField)
    testWindDen =flt.getFLHDenSingle(bField,AWind)
    twistDensity = flt.twistDen(bField,curlField,0.1)
    
    # set up the field line tracer
    
    grid_spacing = [20/192,20/192,20/192]

    #set domain lengths

    lx = grid_spacing[0]*192
    ly = grid_spacing[1]*192
    lz = grid_spacing[2]*192

    #set number of points used to calculate the distribution from
    nx = 200
    ny = 200

    # set a minimum strength of field line cut off 
    bCut = 0.01

    # set z value for anchoring plane (use photosphere z index = 116 here)
    zv = 0
    
    #interpolate the field for tracer
    BxInterp,ByInterp,BzInterp = flt.getInterpolatedFieldSingle(bField,grid_spacing[0],grid_spacing[1],grid_spacing[2])
    
    #calculate and prepare curves
    fieldLinesList,goodSeeds,seeds = flt.prepareCurves(bField,BxInterp,ByInterp,BzInterp,grid_spacing,[6,14],[6,14],zv,nx,ny,bCut)
    
    #interpolate the quantities of interest
    flhInterp = flt.getInterpolatedQuantity(testFLHDen,grid_spacing)
    flhBzInterp = flt.getInterpolatedQuantity(testFLHDen*bz[:,:,0],grid_spacing)
    flwindInterp = flt.getInterpolatedQuantity(testWindDen,grid_spacing)
    twistInterp = flt.getInterpolatedQuantity(twistDensity,grid_spacing)
    
    flh,indexesflh = flt.fieldLineIntegratedQuantity(flhInterp,goodSeeds,fieldLinesList,seeds,nx,ny)
    flhBz,indexesflh = flt.fieldLineIntegratedQuantity(flhBzInterp,goodSeeds,fieldLinesList,seeds,nx,ny)
    flw,indexesflw = flt.fieldLineIntegratedQuantity(flwindInterp,goodSeeds,fieldLinesList,seeds,nx,ny)
    twistF,indexestwist = flt.fieldLineIntegratedQuantity(twistInterp,goodSeeds,fieldLinesList,seeds,nx,ny)

