# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 14:42:56 2022

@author: uqczhan2
if hydraulic conductivity of the sandy layer
is reduced, the large hydraulic conductivity in the oscillating boundary
needs to be redued accordingLy_m 
"""

import os
import sys
import glob
import platform
import shutil
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

# run installed version of flopy or add local path
try:
    import flopy
except:
    fpth = os.path.abspath(os.path.join('..', '..'))
    sys.path.append(fpth)
    import flopy

print(sys.version)
print('numpy version: {}'.format(np.__version__))
print('matplotlib version: {}'.format(mpl.__version__))
print('pandas version: {}'.format(pd.__version__))
print('flopy version: {}'.format(flopy.__version__))

import flopy
import numpy as np
import matplotlib.pyplot as plt



import matplotlib.pylab as pylab
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (15, 5),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}

# %% model set_up
Lx_m   = 300.    # from plot, y is plotted from left to right
Ly_m   = 300.     # from plot, y is plotted upward
ztop_m = 0.    # top elevation of z axis (gravity)
zbot_m = -20.  # bottom elevation of z axis (gravity)
nlay   = 1     # number of layers (gravity)
nrow   = 20     # number of rows
ncol   = 20
delr_m = Lx_m / ncol
delc_m = Ly_m / nrow
delv_m = (ztop_m - zbot_m) / nlay             # layer thickness
vertices_elev_layer_l_list_m = np.linspace(ztop_m, zbot_m, nlay + 1)  # layer elevation array


#vertices_x_coordinate_r_list_m = arrange(0,1000,10)

#dis.itmuni_dict
#dis.itmuni
#hk = 0.3   # how to change the time unit in flopy?  ITMUNI in dis #dis.itmuni_dict horizontal hydraulic conductivity
#vka = 0.3  # vertical hydraulic conductivity
hk_lrc_list        = np.ones((nlay, nrow, ncol), dtype=np.int32) * 10. #*30.   # making hydraulic conductivity array  [lay_row_column]
vka_lrc_list       = hk_lrc_list


sy = 0.25    # specific yield
ss = 1.e-4   #  specific storitivity 
laytyp_l_list = np.ones(nlay)   # vunconfined 1 

#%% the ibound will be defined different this time where we start with 1-D array and later move to 3-D


# interesting to see that the it is better to be converged, if the head is above zero 
strt_lrc_list = 1. * np.ones((nlay, nrow, ncol), dtype=np.float32)   # initial hydraulic head

# time step parameters
nper   = 3                  # number of stress periods
perlen = [100, 100, 100]  # length of stress periods  days 
nstp   = [100, 100, 100]   # number of steps per stress periods
steady = [False, False, False]   # if the results needs to be in steady state.

x_coordiate_r_list_m = np.arange(0,Lx_m +delr_m ,delr_m)
x_coordiate_cell_r_list_m = np.arange(delr_m/2,Lx_m,delr_m)
cell_thickness_r_list = 2 * np.pi * x_coordiate_cell_r_list_m

#%%
modelname = 'ewatering'
mf = flopy.modflow.Modflow(modelname, exe_name='mf2005')
dis = flopy.modflow.ModflowDis(mf, 
                               nlay,                    # number of layers of cells, not layer of vertices
                               nrow, 
                               ncol,            # the is the number of cells, not the number of vertices
                               delr   = delr_m, #np.arange(ncol),#  delr_m, #x_coordiate_cell_r_list_m, #delr_m, #delr (float or array of floats (ncol), optional) – An array of spacings along a row (the default is 1.0).
                               delc   = delc_m, #np.arange(nrow),# delc (float or array of floats (nrow), optional) – An array of spacings along a column (the default is 0.0).
                               top    = ztop_m, # An array of the top elevation of layer 1ztop_m, 
                               botm   = vertices_elev_layer_l_list_m[1:],
                               nper   = nper, 
                               perlen = perlen, 
                               nstp   = nstp, 
                               steady = steady,
                               itmuni = 4)      #itmuni 4 means time unit is days
# https://github.com/modflowpy/flopy/blob/5fcf9709ec8a655106c57d55657c47b1d4987812/examples/Notebooks/flopy3_gridgen.ipynb


    
#%% use of gridgen 
#https://github.com/modflowpy/flopy/blob/5fcf9709ec8a655106c57d55657c47b1d4987812/examples/Notebooks/flopy3_gridgen.ipynb
# setup the active domain
from flopy.utils.gridgen import Gridgen 
# Check and make sure the data folder exists.
model_ws = os.path.join('.', 'data')
if not os.path.exists(model_ws):
    os.makedirs(model_ws)
gridgen_ws = os.path.join(model_ws, 'gridgen')
if not os.path.exists(gridgen_ws):
    os.makedirs(gridgen_ws)
print('Model workspace is : {}'.format(model_ws))
print('Gridgen workspace is : {}'.format(gridgen_ws))


#%%  define a circle in the middle
coordinate_centre_circle_xy_m =  (Lx_m/2.0,Ly_m/2.0)
radius_circle_m   = 120

poly_circle_xy_list=[[[]]]

for i in np.arange(0, 2*np.pi, np.pi/30):
    
    poly_circle_xy_list[0][0].append((coordinate_centre_circle_xy_m[0] + radius_circle_m * np.cos(i),
                                      coordinate_centre_circle_xy_m[1] + radius_circle_m * np.sin(i)
                                      ))
    
# the below line is needed as the fist and final point needs to be exactly (meaning 5~=4.999999999) the same.
poly_circle_xy_list[0][0].append(poly_circle_xy_list[0][0][0])
#%%


g = Gridgen(dis, model_ws=gridgen_ws)

g.build()
#gridgen_ws = os.path.join(model_ws, 'gridgen')

adshp = os.path.join(gridgen_ws, 'ad0')

adpoly = [[[(0, 0), (0, 60), (40, 80), (60, 0), (0, 0)]]]
adline = [[[(0,0),(Lx_m,0),(Lx_m,Ly_m),(0,Ly_m),(0,0)]]]
# g.add_active_domain(adpoly, range(nlay))

adpoly_intersect = g.intersect(poly_circle_xy_list, 'polygon', 0)  # the number at the third argument refers to the layers.
adline_intersect = g.intersect(adline,'line',0)
#adpoly_intersect = g.intersect(poly_circle_xy_list, 'polygon', 1)
print(adpoly_intersect.dtype.names)
print(adpoly_intersect)
print(adpoly_intersect.nodenumber)

print(adline_intersect.nodenumber)


#a = np.zeros((g.nodes), dtype=int)
#a = np.zeros((ncol*nrow),dtype=int) + 2 
ibound_1d_list = np.zeros((ncol*nrow), dtype=int) + 3  # active cell
rf2shp = os.path.join(gridgen_ws, 'rf2')

#rf2shp = os.path.join(gridgen_ws, 'rf2')
#%%  plot IBOUND
#a[adpoly_intersect.nodenumber] = 2
ibound_1d_list[adpoly_intersect.nodenumber]  = 2
ibound_1d_list[adline_intersect.nodenumber]  = -1  # ibound of -1 means the hydraulic head will be stick into the original value.
ibound_lrc_list = ibound_1d_list.reshape(nlay,nrow,ncol)

fig = plt.figure(figsize=(15, 15))
ax = fig.add_subplot(1, 1, 1, aspect='equal')
arr=g.plot(ax, a=ibound_1d_list, masked_values=[0], edgecolor='none', cmap='jet')
mm = flopy.plot.PlotMapView(model=mf)
mm.plot_grid()
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
ax.set_xlabel('X (m)', fontsize=40)
ax.set_ylabel('Y (m)', fontsize=40)
ax.set_title('IBOUND', fontsize=50)
#ax.colorbar(shrink=0.5, ax=ax)
#plt.colorbar(cax=ax)
cbar=plt.colorbar(arr, shrink=0.8, ax=ax, )
cbar.ax.tick_params(labelsize=30)
#quadmesh = mf.plot_ibound() 
#flopy.plot.plot_shapefile(rf2shp, ax=ax, facecolor='yellow', alpha=0.25)


#%% add bas, lpf and pcg package
bas = flopy.modflow.ModflowBas(mf, 
                               ibound = ibound_lrc_list, 
                               strt   = strt_lrc_list
                               )
# https://modflowpy.github.io/flopydoc/mflpf.html
lpf = flopy.modflow.ModflowLpf(mf, 
                               hk     = hk_lrc_list, 
                               vka    = vka_lrc_list, 
                               sy     = sy, 
                               ss     = ss, 
                               laytyp = laytyp_l_list, 
                               ipakcb = 53    ,
                               hdry   = +1e-30,
                               wetfct = 0.1   ,
                               iwetit = 3     ,
                               laywet = 1     ,
                               wetdry = -1    )

pcg = flopy.modflow.ModflowPcg(mf)


#%% plot grid and ibound to show results plot the vertical vview of the model
fig = plt.figure(figsize=(12, 9))
ax  = fig.add_subplot(1, 1, 1, aspect="equal")
modeLx_msect = flopy.plot.PlotCrossSection(model=mf, line={'Row': 10})
arr = modeLx_msect.plot_array(ibound_lrc_list[0,:,:])
modeLx_msect.plot_grid()
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
ax.set_xlabel('X (m)', fontsize=40)
ax.set_ylabel('Y (m)', fontsize=40)
ax.set_title('IBOUND', fontsize=50)
#ax.colorbar(shrink=0.5, ax=ax)
#plt.colorbar(cax=ax)
cbar=plt.colorbar(arr, shrink=0.2, ax=ax)
cbar.ax.tick_params(labelsize=30)

# %%

#dis.sr.xcentergrid
#dis.sr.ycentergrid
#dis.sr.xgrid
#dis.sr.ygrid
#np.ma.masked_equal(dis.sr.ygrid,1000) # very useful command to find specific file locations
#modelmap.sr.vertices
#flopy.plot.plotutil.cell_value_points
#
#
#modeLx_msect = flopy.plot.Modelc_mrossSection(model=mf, line={'Row': 0})
#modeLx_msect.elev
#
#modeLx_msect.dis
#modeLx_msect.xpts
#modeLx_msect.xcentergrid
#modeLx_msect.zcentergrid



#chd={0:[
#       [0,0,0,1,2],
#       [0,0,2,1,2],
#       ]
#    }

ibound_chd_mask=np.ma.masked_equal(ibound,ibound_chd)

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


#flopy.modflow.mfchd.ModflowChd(model=mf,stress_period_data=chd)

stress_period_data = {}
for kper in range(nper):
    for kstp in range(nstp[kper]):
        if np.mod(kstp,49)==0:
            stress_period_data[(kper, kstp)] = ['save head',
                                                'save drawdown',
                                                'save budget',
                                                'print head',
                                                'print budget']
            
oc = flopy.modflow.ModflowOc(mf, stress_period_data=stress_period_data,compact=True)


mf.write_input()


# %%
# Run the model
success, mfoutput = mf.run_model(silent=True, pause=False, report=True)
if not success:
        raise Exception('MODFLOW did not terminate normalLy_m.')


# Imports
import matplotlib.pyplot as plt
import flopy.utils.binaryfile as bf


# Create the headfile and budget file objects
headobj = bf.HeadFile(modelname+'.hds')
times_headobj = headobj.get_times()
cbbobj = bf.CellBudgetFile(modelname+'.cbc')
times_cbbobj = cbbobj.get_times()
# beginning of the first stress period
# end of the first stress period
#times_output_list_day = [100.0, 101.0, 201.0,251]  
#times_output_list_day = [100/2, 200.0/2, 300/2] 
#times_output_list_day = [100, 200.0, 300] 
#times_output_list_day = [50, 150.0, 250] 
times_output_list_day = [99, 101.0, 201,250]

fig = plt.figure(figsize=(8, 3))
ax = fig.add_subplot(1, 1, 1)
#modeLx_msect = flopy.plot.Modelc_mrossSection(model=mf, line={'Column': 5})  # this will onLy_m work when nrow is more than 1
##CM modeLx_msect = flopy.plot.Modelc_mrossSection(model=mf, line={'Row': 0})
modeLx_msect = flopy.plot.PlotCrossSection(model=mf, line={'Row': 0})
patches = modeLx_msect.plot_ibound()
linecollection = modeLx_msect.plot_grid()
t = ax.set_title('Row 0 Cross-Section with IBOUND Boundary Conditions')

#plot(dis.sr.xcentergrid.shape,head[0,0,:])
# head = headobj.get_data(totim=times_output_list_day[2])
# ax.plot(dis.sr.xcentergrid[0,:],head[-1,0,:])
head = headobj.get_data(totim=times_output_list_day[1])
ax.plot(dis.get_node_coordinates()[1],head[-1,0,:])
head = headobj.get_data(totim=times_output_list_day[1])
head[head==1e-30]=np.nan



# %% Plot the Head distribution at the end of stress period 1
#
fig = plt.figure(figsize=(8, 3))
ax  = fig.add_subplot(1, 1, 1)
t   = ax.set_title('Head distribution at the end of stress period 1, day %i' %(times_output_list_day[0]))
head = headobj.get_data(totim=times_output_list_day[0])
modeLx_msect = flopy.plot.PlotCrossSection(model=mf, line={'Row': 0})
arr = modeLx_msect.plot_array(head)
grd = modeLx_msect.plot_grid()
ax.plot(dis.get_node_coordinates()[1] , head[-1,0,:], linewidth=5.0)
plt.colorbar(arr, shrink=1, ax=ax)
ax.set_xlabel('X (m)')
ax.set_ylabel('Z (m)')


times = cbbobj.get_times()
qx = cbbobj.get_data(text="flow right face", totim=times_output_list_day[0])[0]
qy = np.zeros((nlay, nrow, ncol), dtype=float)
qz = cbbobj.get_data(text="flow lower face", totim=times_output_list_day[0])[0]

modeLx_msect.plot_vector(qx, qy, -qz, color="white", kstep=1, hstep=1)
# %% Plot the Head distribution at the end of stress period 1
#
fig = plt.figure(figsize=(8, 3))
ax  = fig.add_subplot(1, 1, 1)
t   = ax.set_title('Head distribution at the end of stress period 1, day %i' %(times_output_list_day[1]))
head = headobj.get_data(totim=times_output_list_day[1])
modeLx_msect = flopy.plot.PlotCrossSection(model=mf, line={'Row': 0})
arr = modeLx_msect.plot_array(head)
grd = modeLx_msect.plot_grid()
ax.plot(dis.get_node_coordinates()[1] , head[-1,0,:], linewidth=5.0)
plt.colorbar(arr, shrink=1, ax=ax)
ax.set_xlabel('X (m)')
ax.set_ylabel('Z (m)')

times = cbbobj.get_times()
qx = cbbobj.get_data(text="flow right face", totim=times_output_list_day[1])[0]
qy = np.zeros((nlay, nrow, ncol), dtype=float)
qz = cbbobj.get_data(text="flow lower face", totim=times_output_list_day[1])[0]

modeLx_msect.plot_vector(qx, qy, -qz, color="white", kstep=1, hstep=1)


# %% Plot the Head distribution at the end of stress period 2
#
fig = plt.figure(figsize=(8, 3))
ax  = fig.add_subplot(1, 1, 1)
t   = ax.set_title('Head distribution at the end of stress period 1, day %i' %(times_output_list_day[2]))
head = headobj.get_data(totim=times_output_list_day[2])
modeLx_msect = flopy.plot.PlotCrossSection(model=mf, line={'Row': 0})
arr = modeLx_msect.plot_array(head)
grd = modeLx_msect.plot_grid()
ax.plot(dis.get_node_coordinates()[1], head[-1,0,:], linewidth=5.0,)
plt.colorbar(arr, shrink=1, ax=ax)
ax.set_xlabel('X (m)')
ax.set_ylabel('Z (m)')

times = cbbobj.get_times()
qx = cbbobj.get_data(text="flow right face", totim=times[4])[0]
qy = np.zeros((nlay, nrow, ncol), dtype=float)
qz = cbbobj.get_data(text="flow lower face", totim=times[4])[0]

modeLx_msect.plot_vector(qx, qy, -qz, color="white", kstep=1, hstep=1)
#CMsurf = ax.plot_surface(modeLx_msect.xcenters, modeLx_msect.xcenters, head[:,0,:], cmap=cm.coolwarm,
#CM                       linewidth=0, antialiased=False)
# %% Plot the Head distribution at the end of stress period 3
#
fig = plt.figure(figsize=(8, 3))
ax  = fig.add_subplot(1, 1, 1)
t   = ax.set_title('Head distribution at the end of stress period 1, day %i' %(times_output_list_day[3]))
head = headobj.get_data(totim=times_output_list_day[3])
modeLx_msect = flopy.plot.PlotCrossSection(model=mf, line={'Row': 0})
arr = modeLx_msect.plot_array(head)
grd = modeLx_msect.plot_grid()
ax.plot(dis.get_node_coordinates()[1] , head[-1,0,:], linewidth=5.0,)
plt.colorbar(arr, shrink=1, ax=ax)
ax.set_xlabel('X (m)')
ax.set_ylabel('Z (m)')

times = cbbobj.get_times()
qx = cbbobj.get_data(text="flow right face", totim=times[5])[0]
qy = np.zeros((nlay, nrow, ncol), dtype=float)
qz = cbbobj.get_data(text="flow lower face", totim=times[5])[0]

modeLx_msect.plot_vector(qx, qy, -qz, color="white", kstep=1, hstep=1)
#CMsurf = ax.plot_surface(modeLx_msect.xcenters, modeLx_msect.xcenters, head[:,0,:], cmap=cm.coolwarm,
#CM                       linewidth=0, antialiased=False)
