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

# print the versions used to run this model, a important information to
# reproduce the result a couple of years later
print(sys.version)
print('numpy version: {}'.format(np.__version__))
print('matplotlib version: {}'.format(mpl.__version__))
print('pandas version: {}'.format(pd.__version__))
print('flopy version: {}'.format(flopy.__version__))

# parameters for plotting
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
# time step parameters
nper   = 2                  # number of stress periods
perlen = [100,100]          # length of stress periods  days 
nstp   = [100,100]          # number of steps per stress periods
tsmult = 1. # the multiplier for the length of successive time steps, 1 means time step is always the same
bool_steady_state_stress_period = [False, False]     # if the results needs to be in steady state
step_interval_output = 2    # output will be saved every # of intervals

#vertices_x_coordinate_r_list_m = arrange(0,1000,10)

#dis.itmuni_dict
#dis.itmuni
#hk = 0.3   # how to change the time unit in flopy?  ITMUNI in dis #dis.itmuni_dict horizontal hydraulic conductivity
#vka = 0.3  # vertical hydraulic conductivity
hk_lrc_list        = np.ones((nlay, nrow, ncol), dtype=np.int32) * 10. #*30.   # making hydraulic conductivity array  [lay_row_column]
vka_lrc_list       = hk_lrc_list


sy = 0.25    # specific yield, equivalent to effective porosity
ss = 1.e-4   #  specific storitivity, corresponding to the compressibity of solid matrix and fluids (water)
laytyp_l_list = np.ones(nlay)   # vunconfined 1 



#%%
modelname = 'ewatering'
#mf = flopy.modflow.Modflow(modelname, exe_name='mf2005')
ws = os.path.abspath('../Model/')
sim = flopy.mf6.MFSimulation(sim_name=modelname, 
                             version='mf6',
#                             exe_name='../Exe/mf6',  #comment this line means flopy will look for mf6 from system folders
                             sim_ws=ws)

tdis = flopy.mf6.ModflowTdis(sim, 
                             time_units='DAYS',
                             nper=nper, 
                             perioddata=[[perlen[0],nstp[0],tsmult],[perlen[1],nstp[1],tsmult]]
                            )

#%% dis
fModName = 'FlowModel'
gwf = flopy.mf6.ModflowGwf(sim, 
                           modelname  =  fModName, 
                           newtonoptions  = True)
#%%  ims package
nouter, ninner = 700, 300
hclose, rclose, relax = 1e-8, 1e-6, 0.97   # convergene criteria
imsgwf = flopy.mf6.ModflowIms(sim, print_option='ALL',
                                  outer_dvclose=hclose,
                                  outer_maximum=nouter,
                                  under_relaxation='NONE',
                                  inner_maximum=ninner,
                                  inner_dvclose=hclose, 
                                  rcloserecord=rclose,
                                  linear_acceleration='BICGSTAB',
                                  scaling_method='NONE',
                                  reordering_method='NONE',
                                  relaxation_factor=relax,
                                  filename='{}.ims'.format(fModName))

#%%

idomain = np.full((nlay, nrow, ncol), 1) # similar to IBOUND where 0 means inactive cell, 1 means active cell
# idomain[0, 0, 1:6] = 0
# idomain[1, 0, 2:5] = 0
# idomain[2, 0, 3:4] = 0
dis = flopy.mf6.ModflowGwfdis(gwf, 
                              nlay=nlay, 
                              nrow=nrow, 
                              ncol=ncol,
                              delr=delr_m, 
                              delc=delc_m,
                              top= ztop_m, 
                              botm=vertices_elev_layer_l_list_m[1:], 
                              idomain=1)



# dis = flopy.modflow.ModflowDis(mf, 
#                                nlay,                    # number of layers of cells, not layer of vertices
#                                nrow, 
#                                ncol,            # the is the number of cells, not the number of vertices
#                                delr   = delr_m, #np.arange(ncol),#  delr_m, #x_coordiate_cell_r_list_m, #delr_m, #delr (float or array of floats (ncol), optional) – An array of spacings along a row (the default is 1.0).
#                                delc   = delc_m, #np.arange(nrow),# delc (float or array of floats (nrow), optional) – An array of spacings along a column (the default is 0.0).
#                                top    = ztop_m, # An array of the top elevation of layer 1ztop_m, 
#                                botm   = vertices_elev_layer_l_list_m[1:], # An array of the bottom elevation of layer (cell)  
#                                nper   = nper, 
#                                perlen = perlen, 
#                                nstp   = nstp, 
#                                steady = bool_steady_state_stress_period,
#                                itmuni = 4)      # itmuni 4 means time unit is days for this model
# # https://github.com/modflowpy/flopy/blob/5fcf9709ec8a655106c57d55657c47b1d4987812/examples/Notebooks/flopy3_gridgen.ipynb


    
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
    
    poly_circle_xy_list[0][0].append((coordinate_centre_circle_xy_m[0] + 
                                      radius_circle_m * np.cos(i),
                                      coordinate_centre_circle_xy_m[1] + 
                                      radius_circle_m * np.sin(i)
                                      ))
    
# the below line is needed as the fist and final point needs to be exactly (meaning 5~=4.999999999) the same.
poly_circle_xy_list[0][0].append(poly_circle_xy_list[0][0][0])
#%%
g = Gridgen(dis, model_ws=gridgen_ws)
g.build()

#gridgen_ws = os.path.join(model_ws, 'gridgen')

adshp = os.path.join(gridgen_ws, 'ad0')

#adpoly = [[[(0, 0), (0, 60), (40, 80), (60, 0), (0, 0)]]]
adline = [[[(0,0),(Lx_m,0),(Lx_m,Ly_m),(0,Ly_m),(0,0)]]]
# g.add_active_domain(adpoly, range(nlay))

adpoly_intersect = g.intersect(poly_circle_xy_list, 'polygon', 0)  # the number at the third argument refers to the layers.
adline_intersect = g.intersect(adline,'line',0)
#adpoly_intersect = g.intersect(poly_circle_xy_list, 'polygon', 1)
print(adpoly_intersect.dtype.names)
print(adpoly_intersect)
print(adpoly_intersect.nodenumber)
print(adline_intersect.nodenumber)

#g.add_refinement_features(poly_circle_xy_list, 'polygon', 1, range(nlay))

#a = np.zeros((g.nodes), dtype=int)
#a = np.zeros((ncol*nrow),dtype=int) + 2 
ibound_1d_list = np.zeros((ncol*nrow), dtype=int) + 3  # active cell
rf2shp = os.path.join(gridgen_ws, 'rf0')
#rf2shp = os.path.join(gridgen_ws, 'rf2')
#%%  plot IBOUND
#a[adpoly_intersect.nodenumber] = 2
ibound_chd = 2            # all the ibounds that will be subjected to recharge will be tagged with 2.
ibound_constant_head = -1 # ibound of -1 means the hydraulic head will be stick into the original value.
ibound_1d_list[adpoly_intersect.nodenumber]  = ibound_chd
ibound_1d_list[adline_intersect.nodenumber]  = ibound_constant_head
ibound_lrc_list = ibound_1d_list.reshape(nlay,nrow,ncol)

dis.idomain = ibound_lrc_list

fig = plt.figure(figsize=(15, 15))
ax = fig.add_subplot(1, 1, 1, aspect='equal')
arr=g.plot(ax, a=dis.idomain, masked_values=[0], edgecolor='none', cmap='jet')
mm = flopy.plot.PlotMapView(model=gwf)
#flopy.plot.plot_shapefile(rf2shp, ax=ax, facecolor='yellow', alpha=0.25)
mm.plot_grid()
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
ax.set_xlabel('X (m)', fontsize=40)
ax.set_ylabel('Y (m)', fontsize=40)
ax.set_title('IDOMAIN', fontsize=50)
#ax.colorbar(shrink=0.5, ax=ax)
#plt.colorbar(cax=ax)
#cbar=plt.colorbar(arr, shrink=0.8, ax=ax)
#cbar.ax.tick_params(labelsize=30)
#quadmesh = mf.plot_ibound() 
#flopy.plot.plot_shapefile(rf2shp, ax=ax, facecolor='yellow', alpha=0.25)


#%% the ibound will be defined different this time where we start with 1-D array and later move to 3-D


# interesting to see that the it is better to be converged, if the head is above zero 
strt_lrc_list_m = -5. * np.ones((nlay, nrow, ncol), dtype=np.float32)   # initial hydraulic head
ic = flopy.mf6.ModflowGwfic(gwf, strt = strt_lrc_list_m)

# bas = flopy.modflow.ModflowBas(mf, 
#                                ibound = ibound_lrc_list, 
#                                strt   = strt_lrc_list_m
#                                )
#%% #k33 ([double]) –
#k33 (double) is the hydraulic conductivity of the third ellipsoid axis (or the ratio of K33/K if the K33OVERK option is specified); for an unrotated case, this is the vertical hydraulic conductivity.
#When anisotropy is applied, K33 corresponds to the K33 tensor component. All included cells (IDOMAIN > 0) must have a K33 value greater than zero.
# icelltype (integer) flag for each cell that specifies how saturated thickness is treated. 0 means saturated thickness is held constant; 
#:math:`>`0 means saturated thickness varies with computed head when head is below the cell top; 
#:math:`<`0 means saturated thickness varies with computed head unless the THICKSTRT option is in effect. 
# When THICKSTRT is in effect, a negative value of icelltype indicates that saturated thickness will be computed as STRT-BOT and held constant.
# Kh = 1.
# Kv = 1.
npf = flopy.mf6.ModflowGwfnpf(gwf, xt3doptions = False,
                                  save_flows = True,
                                  save_specific_discharge = True,
                                  icelltype = 1,   # meaning that transmissivity changes with heads
                                  k   = hk_lrc_list, 
                                  k33 = vka_lrc_list)

sto = flopy.mf6.ModflowGwfsto(gwf, 
                              sy = sy, 
                              ss = ss, 
                              iconvert=1
                              )
# # https://modflowpy.github.io/flopydoc/mflpf.html
# lpf = flopy.modflow.ModflowLpf(mf, 
#                                hk     = hk_lrc_list, 
#                                vka    = vka_lrc_list, 
#                                sy     = sy, 
#                                ss     = ss, 
#                                laytyp = laytyp_l_list, 
#                                ipakcb = 53    ,
#                                hdry   = +1e-30,
#                                wetfct = 0.1   ,
#                                iwetit = 3     ,
#                                laywet = 1     ,
#                                wetdry = -1    )
#
#%%
# CZ220303 removed 'steps' to make it work
oc = flopy.mf6.ModflowGwfoc(gwf,
                            budget_filerecord ='{}.cbc'.format(fModName),
                            head_filerecord   ='{}.hds'.format(fModName),
                            headprintrecord   =[
                                ('COLUMNS', 10, 'WIDTH', 15,
                                 'DIGITS', 6, 'GENERAL')],
                            saverecord=[('HEAD', 'ALL'),
                                        ('BUDGET', 'ALL')],
                            printrecord=[('HEAD', 'LAST'),
                                         ('BUDGET', 'LAST')])
#pcg = flopy.modflow.ModflowPcg(mf)



#%% plot grid and ibound to show results plot the vertical vview of the model
fig = plt.figure(figsize=(12, 9))
ax  = fig.add_subplot(1, 1, 1, aspect="equal")
modeLx_msect = flopy.plot.PlotCrossSection(model=gwf, line={'Row': 10})
arr = modeLx_msect.plot_array(dis.idomain[0,:,:])
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


# stress_period_data = {}
# for kper in range(nper):
#     for kstp in range(nstp[kper]):
#         if np.mod(kstp,step_interval_output) == 0:
#             stress_period_data[(kper, kstp)] = ['save head',
#                                                 'save drawdown',
#                                                 'save budget',
#                                                 'print head',
#                                                 'print budget']
# # output control            
# oc = flopy.modflow.ModflowOc(mf, 
#                              stress_period_data=stress_period_data,
#                              compact=True)
# recharge package
recharge_rate_mPday = {"0":0.1,"1":0.1}   #0.1 #0.01

# nrchop (int) – is the recharge option code. 
# 1: Recharge to top grid layer only 
# 2: Recharge to layer defined in irch 
# 3: Recharge to highest active cell (default is 3).
nrchop = 1


# rch=flopy.modflow.ModflowRch(mf,
#                              rech   = recharge_rate_mPday,
#                              nrchop = nrchop,
#                              ipakcb = 1,
#                              stress_period_data   = {'0': adpoly_intersect.nodenumber , 
#                                        '1': adpoly_intersect.nodenumber}
#                              )
rch = flopy.mf6.ModflowGwfrcha(
    gwf, 
    recharge   = 0.1,
    #nrchop = 1
    )
#%%
#gwf.write_input()
sim.write_simulation(silent=False)

# %% Run the model
try:
    os.remove(os.path.join(model_ws, "{0}.hds".format(modelname)))
except:
    pass

success,mfoutput = sim.run_simulation(silent=False, pause=False,report=True) #mfoutput = mf.run_model(silent=True, pause=False, report=True)

if not success:
        raise Exception('MODFLOW did not terminate normalLy_m.')

# %% Extracting data
import flopy.utils.binaryfile as bf


# Create the headfile and budget file objects
headobj       = bf.HeadFile(modelname+'.hds')
times_headobj = headobj.get_times()
cbbobj        = bf.CellBudgetFile(modelname+'.cbc')
times_cbbobj  = cbbobj.get_times()
# beginning of the first stress period



# %%
times_output_list_day = [11, 101.0, 201, 301]

fig = plt.figure(figsize=(8, 3))
ax = fig.add_subplot(1, 1, 1)
#modeLx_msect = flopy.plot.Modelc_mrossSection(model=mf, line={'Column': 5})  
# this will onLy_m work when nrow is more than 1
##CM modeLx_msect = flopy.plot.Modelc_mrossSection(model=mf, line={'Row': 0})
modeLx_msect = flopy.plot.PlotCrossSection(model=gwf, line={'Row': int(nrow/2) })
patches = modeLx_msect.plot_ibound()
linecollection = modeLx_msect.plot_grid()
t = ax.set_title('Row 0 Cross-Section with IBOUND Boundary Conditions')

#plot(dis.sr.xcentergrid.shape,head[0,0,:])
# head = headobj.get_data(totim=times_output_list_day[2])
# ax.plot(dis.sr.xcentergrid[0,:],head[-1,0,:])
head = headobj.get_data(totim=times_output_list_day[1])
ax.plot(gwf.modelgrid.xycenters[0],head[-1,int(nrow/2),:])   
head = headobj.get_data(totim=times_output_list_day[1])
head[head==1e-30] = np.nan



# %% Plot the Head distribution at the end of stress period 1
#
fig = plt.figure(figsize=(8, 3))
ax  = fig.add_subplot(1, 1, 1)
t   = ax.set_title('Head distribution at the end of stress period 1, day %i' %(times_output_list_day[0]))
head = headobj.get_data(totim=times_output_list_day[0])
modeLx_msect = flopy.plot.PlotCrossSection(model=gwf, line={'Row': 0})
arr = modeLx_msect.plot_array(head)
grd = modeLx_msect.plot_grid()
ax.plot(gwf.modelgrid.xycenters[0] , head[-1,0,:], linewidth=5.0)
plt.colorbar(arr, shrink=1, ax=ax)
ax.set_xlabel('X (m)')
ax.set_ylabel('Z (m)')


# times = cbbobj.get_times()
# qx = cbbobj.get_data(text="flow right face", totim=times_output_list_day[0])[0]
# qy = np.zeros((nlay, nrow, ncol), dtype=float)
# qz = cbbobj.get_data(text="flow lower face", totim=times_output_list_day[0])[0]

# modeLx_msect.plot_vector(qx, qy, -qz, color="white", kstep=1, hstep=1)
# %% Plot the Head distribution at the end of stress period 1
#
fig = plt.figure(figsize=(8, 3))
ax  = fig.add_subplot(1, 1, 1)
t   = ax.set_title('Head distribution at the end of stress period 1, day %i' %(times_output_list_day[1]))
head = headobj.get_data(totim=times_output_list_day[1])
modeLx_msect = flopy.plot.PlotCrossSection(model=gwf, line={'Row': 0})
arr = modeLx_msect.plot_array(head)
grd = modeLx_msect.plot_grid()
ax.plot(gwf.modelgrid.xycenters[0] , head[-1,0,:], linewidth=5.0)
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
