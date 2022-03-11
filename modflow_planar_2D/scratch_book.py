# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 14:51:37 2022

@author: uqczhan2
"""

#%%  extra script to use
import flopy.discretization as fgrid
import shapely

xoff = 0.
yoff = 0.
angrot = 0.
delc = 10 * np.ones(10, dtype=float)
delr = 10 * np.ones(10, dtype=float)
sgr = fgrid.StructuredGrid(delc, delr, top=None, botm=None, xoff=xoff, yoff=yoff, angrot=angrot)

sgr.plot()

from shapely.geometry import Polygon
p = Polygon(shell=[(15, 15), (20, 50), (35, 80.), (80, 50), 
                   (80, 40), (40, 5), (15, 12)])
ix = GridIntersect(sgr, method="vertex", rtree=True)

#%timeit ix.intersect(p)




#%%
# First step is to set up the plot
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(1, 1, 1, aspect="equal")

# Next we create an instance of the PlotMapView class
mapview = flopy.plot.PlotMapView(model=mf)

# Then we can use the plot_grid() method to draw the grid
# The return value for this function is a matplotlib LineCollection object,
# which could be manipulated (or used) later if necessary.
linecollection = mapview.plot_grid()

t = ax.set_title("Model Grid")



'''
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
fig = plt.figure()
ax = fig.gca(projection='3d')

# Make data.
X = np.arange(-5, 5, 0.25)
Y = np.arange(-5, 5, 0.25)
X, Y = np.meshgrid(X, Y)
R = np.sqrt(X**2 + Y**2)
Z = np.sin(R)
# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)
# Customize the z axis.
ax.set_zlim(-1.01, 1.01)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()
'''

# we are supposed to use CHD boundary here not GHB boundary
'''
# Make list for stress period 1
stageleft = 10.
stageright = 10.
bound_sp1 = []
for il in range(nlay):
    condleft = hk * (stageleft - zbot_m) * delc_m
    condright = hk * (stageright - zbot_m) * delc_m
    for ir in range(nrow):
        bound_sp1.append([il, ir, 0, stageleft, condleft])
        bound_sp1.append([il, ir, ncol - 1, stageright, condright])
print('Adding ', len(bound_sp1), 'GHBs for stress period 1.')

# Make list for stress period 2
stageleft = 10.
stageright = 0.
condleft = hk * (stageleft - zbot_m) * delc_m
condright = hk * (stageright - zbot_m) * delc_m
bound_sp2 = []
for il in range(nlay):
    for ir in range(nrow):
        bound_sp2.append([il, ir, 0, stageleft, condleft])
        bound_sp2.append([il, ir, ncol - 1, stageright, condright])
print('Adding ', len(bound_sp2), 'GHBs for stress period 2.')

# We do not need to add a dictionary entry for stress period 3.
# Flopy will automaticalLy_m take the list from stress period 2 and appLy_m it
# to the end of the simulation, if necessary
stress_period_data = {0: bound_sp1, 1: bound_sp2}

# Create the flopy ghb object
ghb = flopy.modflow.ModflowGhb(mf, stress_period_data=stress_period_data)
'''

#%%
ibound_chd_mask=np.ma.masked_equal(ibound,ibound_chd)


ibound_chd_mask=np.ma.masked_equal(ibound_lrc,ibound_chd)


chd_node_index=np.where(ibound_chd_mask.mask)


stress_period_data = {}

bound_sp0 = []
for i in np.arange(np.sum(ibound_chd_mask.mask)):
   bound_sp0.append([chd_node_index[0][i],chd_node_index[1][i],chd_node_index[2][i],0,0]  )
bound_sp1=[]
for i in np.arange(np.sum(ibound_chd_mask.mask)):
   bound_sp1.append([chd_node_index[0][i],chd_node_index[1][i],chd_node_index[2][i],-10,-10]  )

bound_sp2=[]
for i in np.arange(np.sum(ibound_chd_mask.mask)):
   bound_sp2.append([chd_node_index[0][i],chd_node_index[1][i],chd_node_index[2][i],0,0]  )



stress_period_data={0:bound_sp0,1:bound_sp1,2:bound_sp2}


# write to chd package
chd=flopy.modflow.mfchd.ModflowChd(model=mf,stress_period_data=stress_period_data)


#%% a very quick way to show the surfae image.
# sourced from https://github.com/modflowpy/flopy/blob/develop/examples/groundwater_paper/Notebooks/uspb.ipynb
c = plt.imshow(grndElv, cmap='jet')
plt.colorbar(c);





# %% extracting data from key locations
# https://github.com/connorcleary/code/blob/9b4af7abfe097e03afffaa07af24d7744a080285/flopy/jupyter/flopy.ipynb

# def get_lrc_from_coordinates(x,y,z,gwf=gwf) :
#     """
#     function to get the lrc of the coordinates. 
#     """
#     from flopy.utils.gridintersect import GridIntersect
#     ix = GridIntersect(gwf.modelgrid, method='vertex')
#     ix.intersects(
#     [row,column]= dis.get_rc_from_node_coordinates(x,y)
#     #layer = dis.get_layer(row,column,z)
#     layer = 0
#     return (layer, row , column)


# # def extract_head_from_xyz(x,y,z,gwf=gwf):
# #     result = {'x':x,'y':y,'z':z}
# #     [result['r'] ,result['c']] =
# #         gwf.modelgrid.get_coords(result['x'] ,result['y'])
    
    
# point_1 = {'x':500,'y':5,'z':-15}
# (point_1['l'],point_1['r'],point_1['c'] )= get_lrc_from_coordinates(
#     point_1['x'],
#     point_1['y'],
#     point_1['z'])
# point_1['head_time_ay']=np.array(
#     [ i[point_1['l'],point_1['r'],point_1['c']] for i in head_array_timeid_lrc_m ])
# point_1 ['label']= 'x={:1.1e}'.format(point_1['x']) +  ', z={:1.1e}'.format(point_1['z'])    