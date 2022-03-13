# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 14:42:56 2022

@author: uqczhan2
if hydraulic conductivity of the sandy layer
is reduced, the large hydraulic conductivity in the oscillating boundary
needs to be redued accordingLy_m 

need to install: 
pip install -U csaps

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
model_ws = os.path.join('.', 'data')


print('Model workspace is : {}'.format(os.getcwd()))

#%% unit conversion
mPmm = 0.001 # convert mm to m


#%% model set_up
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

# # time step parameters
# nper   = 2                  # number of stress periods
# perlen = [100,100]          # length of stress periods  days 

# 
# nstp   = [100,100]          # number of steps per stress periods
# tsmult = [1. , 1.] # the multiplier for the length of successive time steps, 1 means time step is always the same
# bool_steady_state_stress_period = [False, False]     # if the results needs to be in steady state

#last for 200 days, using 200 stress period, each stress period has one day
nper = 200
perlen = np.ones(nper)
nstp   = np.ones(nper,dtype=int)
tsmult = np.ones(nper)
bool_steady_state_stress_period = np.zeros(nper,dtype=bool)



stress_period_end_time_days_ay = np.cumsum(perlen)
 # output will be saved every # of intervals, unused at the moment
step_interval_output = 2   

# period data array input for dis package
dis_perioddata_ay = []
for i in np.arange(nper):
    dis_perioddata_ay.append([perlen[i],nstp[i],tsmult[i]])

hk_mPday = 10.0 
hk_lrc_list        = np.ones((nlay, nrow, ncol), dtype=np.int32) * hk_mPday #*30.   # making hydraulic conductivity array  [lay_row_column]
vka_lrc_list       = hk_lrc_list


sy = 0.25    # specific yield, equivalent to effective porosity
ss = 1.e-4   #  specific storitivity, corresponding to the compressibity of solid matrix and fluids (water)

#%%  read field data, which forms input and are for verification of the model.
field_data_ws = os.path.join(model_ws, 'field_data.csv')
field_data_df = pd.read_csv(field_data_ws,
                            parse_dates=['date_time'],
                            index_col=0,
                            )
# the row index is time rather than 0,1,2. this will help analysis such as 
#getting daily averge, max, min or data at a specific hours. 
# to get data in row 2 use field_data_df.loc(2) 
#field_data_df= field_data_df.set_index('date_time')
field_data_df['time_elapsed_days'] = (field_data_df.index - field_data_df.index[0])/np.timedelta64(1, 'D')


#%% create a time array, starting from zero, the list of times
# this is needed as we need to use surface water level to calcuate the 
# recharge at every step.
# this does not consider the impact of tmulti
# CZ220311 if tdis has this array, we do not need to produce here anymore.
# 
# for i in np.arange(nper):
#     time_ay_current_stress_period = np.linspace(0,perlen[i]-perlen[i]/nstp[i],nstp[i]) 
#     if i == 0: 
#         times_ay_days= time_ay_current_stress_period
#     else:
#         times_ay_days = np.concatenate((times_ay_days, 
#                                         times_ay_days[-1] + 
#                                         perlen[i-1]/nstp[i-1] + 
#                                         time_ay_current_stress_period
#                                         ))
# replaced by gwf.modeltime.totim


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
                             perioddata=dis_perioddata_ay
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
#%% cubic spline interpolating the surface water depth folloing the measurement

from csaps import csaps
csaps_coef = 0.85
surface_water_depth_rch_input_totim_ay_m = csaps(field_data_df['time_elapsed_days'], 
           field_data_df['surface_water_depth_m'],
           gwf.modeltime.totim, 
           smooth = csaps_coef)

#%%
fig=plt.figure()

ax=field_data_df.plot(x = 'time_elapsed_days',
                   y = 'surface_water_depth_m')
plt.plot(gwf.modeltime.totim,
         surface_water_depth_rch_input_totim_ay_m,
         'o')
ax.set_title('interpolate surface water depth for recharge input, coef = ' + str(csaps_coef) )
#%%

idomain = np.full((nlay, nrow, ncol), 1) # similar to idomain where 0 means inactive cell, 1 means active cell
# idomain[0, 0, 1:6] = 0
# idomain[1, 0, 2:5] = 0
# idomain[2, 0, 3:4] = 0
dis = flopy.mf6.ModflowGwfdis(gwf, 
                              nlay = nlay, 
                              nrow = nrow, 
                              ncol = ncol,
                              delr = delr_m, 
                              delc = delc_m,
                              top  = ztop_m, 
                              botm = vertices_elev_layer_l_list_m[1:], 
                              idomain=1   # temperorily presscribed
                              )
# # https://github.com/modflowpy/flopy/blob/5fcf9709ec8a655106c57d55657c47b1d4987812/examples/Notebooks/flopy3_gridgen.ipynb


    
#%% use of gridgen 
#https://github.com/modflowpy/flopy/blob/5fcf9709ec8a655106c57d55657c47b1d4987812/examples/Notebooks/flopy3_gridgen.ipynb
# setup the active domain
from flopy.utils.gridgen import Gridgen 
# Check and make sure the data folder exists.

if not os.path.exists(model_ws):
    os.makedirs(model_ws)
gridgen_ws = os.path.join(model_ws, 'gridgen')
if not os.path.exists(gridgen_ws):
    os.makedirs(gridgen_ws)
print('Model workspace is : {}'.format(model_ws))
print('Gridgen workspace is : {}'.format(gridgen_ws))


#%%  define a circle in the middle
coordinate_centre_circle_xy_m =  (Lx_m/2.0,Ly_m/2.0)
radius_circle_m   = 100

poly_circle_xy_list=[[[]]]

for i in np.arange(0, 2*np.pi, np.pi/30):
    
    poly_circle_xy_list[0][0].append((coordinate_centre_circle_xy_m[0] + 
                                      radius_circle_m * np.cos(i),
                                      coordinate_centre_circle_xy_m[1] + 
                                      radius_circle_m * np.sin(i)
                                      ))
    
# the below line is needed as the fist and final point needs to be exactly (meaning 5~=4.999999999) the same.
poly_circle_xy_list[0][0].append(poly_circle_xy_list[0][0][0])

g = Gridgen(dis, model_ws=gridgen_ws)
g.build()

#gridgen_ws = os.path.join(model_ws, 'gridgen')

adshp = os.path.join(gridgen_ws, 'ad0')


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


idomain_1d_list = np.zeros((ncol*nrow), dtype=int) + 3  # active cell
rf2shp = os.path.join(gridgen_ws, 'rf0')
#rf2shp = os.path.join(gridgen_ws, 'rf2')
#%%  plot idomain
#a[adpoly_intersect.nodenumber] = 2
idomain_rch = 2  # all the idomains that will be subjected to recharge will be tagged with 2.
idomain_constant_head = 5 # chd
idomain_1d_list[adpoly_intersect.nodenumber]  = idomain_rch
idomain_1d_list[adline_intersect.nodenumber]  = idomain_constant_head
idomain_lrc_list = idomain_1d_list.reshape(nlay,nrow,ncol)
dis.idomain = idomain_lrc_list


#ax.colorbar(shrink=0.5, ax=ax)
#plt.colorbar(cax=ax)
#cbar=plt.colorbar(arr, shrink=0.8, ax=ax)
#cbar.ax.tick_params(labelsize=30)
#quadmesh = mf.plot_idomain() 
#flopy.plot.plot_shapefile(rf2shp, ax=ax, facecolor='yellow', alpha=0.25)


#%% the idomain will be defined different this time where we start with 1-D array and later move to 3-D


# interesting to see that the it is better to be converged, if the head is above zero 
strt_m =-5 
strt_lrc_list_m = strt_m * np.ones((nlay, nrow, ncol), dtype=np.float32)   # initial hydraulic head
ic = flopy.mf6.ModflowGwfic(gwf, strt = strt_lrc_list_m)


#%% #k33 ([double]) –
#k33 (double) is the hydraulic conductivity of the third ellipsoid axis (or the ratio of K33/K if the K33OVERK option is specified); for an unrotated case, this is the vertical hydraulic conductivity.
#When anisotropy is applied, K33 corresponds to the K33 tensor component. All included cells (IDOMAIN > 0) must have a K33 value greater than zero.
# icelltype (integer) flag for each cell that specifies how saturated thickness is treated. 0 means saturated thickness is held constant; 
#:math:`>`0 means saturated thickness varies with computed head when head is below the cell top; 
#:math:`<`0 means saturated thickness varies with computed head unless the THICKSTRT option is in effect. 
# When THICKSTRT is in effect, a negative value of icelltype indicates that saturated thickness will be computed as STRT-BOT and held constant.
# k22 is equal to k by default
npf = flopy.mf6.ModflowGwfnpf(gwf,
                              xt3doptions = False,
                              save_flows = True,
                              save_specific_discharge = True,
                              icelltype = 1,   # meaning that transmissivity changes with heads
                              k   = hk_lrc_list, 
                              k33 = vka_lrc_list)

sto = flopy.mf6.ModflowGwfsto(gwf, 
                              sy = sy, 
                              ss = ss, 
                              iconvert = 1
                              )
# # https://modflowpy.github.io/flopydoc/mflpf.html
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


#%% recharge package

#https://stackoverflow.com/questions/7961363/removing-duplicates-in-lists
#list(set(t))
# duplicate nodes exists in the intersect nodenumber, which needs to be removed

rch_spd_dict = {}

recharge_rate_spd_ay = np.zeros(nper,dtype=float)
# during the first 100 days, the recharge rate is 0.1 m/day
recharge_rate_spd_ay [ stress_period_end_time_days_ay <= 100 ] = 0.1

# list of variables used for calculating recharge at each cells
# surface_water_depth_rch_input_totim_ay_m
# gwf.modeltime.totim
# gwf.modelgrid.xcellcenters
# gwf.modelgrid.ycellcenters
# surface_elevation_cell_rl_ay_m




for per in np.arange(nper):
    rch_spd = []
    for i in list(set(adpoly_intersect.nodenumber)) :
        coord= gwf.modelgrid.get_lrc(i)
        rch_spd.append([0, coord[0][1], coord[0][2], recharge_rate_spd_ay[per]] )
    rch_spd_dict[per] =  rch_spd


rch = flopy.mf6.ModflowGwfrch(
    gwf, 
    stress_period_data = rch_spd_dict
)



#%% constant head package
chd_spd = []

boundary_constant_head_m= -5
for i in list(set(adline_intersect.nodenumber)):
    coord= gwf.modelgrid.get_lrc(i)
    chd_spd.append(  [ (0 , coord[0][1], coord[0][2]) , 
                      boundary_constant_head_m  ] )

chd = flopy.mf6.ModflowGwfchd(
    gwf,
    stress_period_data={0: chd_spd,1: chd_spd},
)
    


# %% observation package
# installing wells 
# Note: the row starts from top to bottom, so the first row has 
# a coordinate of 300, the last row has a corrdinate of 0

from flopy.utils.gridintersect import GridIntersect
from shapely.geometry import Point
ix = GridIntersect(gwf.modelgrid, method='vertex') # CZ220310 if the point is at the grid, both cells will be included.

point_centre_pond = Point(150.1,150.1)
point_centre_pond.id = ix.intersects(point_centre_pond)
point_centre_pond.lrc_loc = (0, 
                             point_centre_pond.id['cellids'][0][0], 
                             point_centre_pond.id['cellids'][0][1])  # (layerid, rowid, columnid)

# importing obs id
#obs_wells = r".\objects\obs_wells.csv"
#bs_vertices = np.genfromtxt(obs_wells, skip_header = 1, delimiter = ',')

# intersect with grid
#obs_id = ix.intersects(Point(obs_vertices[0], obs_vertices[1]))
#obs_id

point_SA3 = Point (Point(240.1,150.1) )
point_SA3.cell_id = ix.intersects(point_SA3)
point_SA3.lrc_loc = (0, 
                     point_SA3.cell_id['cellids'][0][0], 
                     point_SA3.cell_id['cellids'][0][1]) # (layerid, rowid, columnid)

point_SA4 = Point (Point(280.1,150.1) )
point_SA4.cell_id = ix.intersects(point_SA4)
point_SA4.lrc_loc = (0, 
                     point_SA4.cell_id['cellids'][0][0], 
                     point_SA4.cell_id['cellids'][0][1]) # (layerid, rowid, columnid)

# dist_to_pond_centre_xy_ay_m = np.zeros([nrow,ncol],dtype=float)

# for i in np.arange(nrow):
#     for j in np.arange(ncol):
#         dist = (gwf.modelgrid.xcellcenters[i][j] - point_centre_pond.x) **
# pond bed slope in tangent value.
pond_bed_slope_tan = 0.003  # np.tan(.5*np.pi/180) np.arctan(0.003)/np.pi*180

dist_to_pond_centre_cell_rl_ay_m = (( gwf.modelgrid.xcellcenters - point_centre_pond.x ) ** 2. + \
                                   ( gwf.modelgrid.ycellcenters - point_centre_pond.y ) ** 2. ) \
                                   ** 0.5

max_depth_m = radius_circle_m * pond_bed_slope_tan

depth_cell_rl_ay_m = max_depth_m - dist_to_pond_centre_cell_rl_ay_m * pond_bed_slope_tan

depth_cell_rl_ay_m [depth_cell_rl_ay_m <  0] =0

surface_elevation_cell_rl_ay_m = ztop_m - depth_cell_rl_ay_m



from mpl_toolkits.mplot3d.axes3d import Axes3D, get_test_data
from matplotlib import cm
#get_ipython().run_line_magic('matplotlib', 'auto')  #allow the graph to pop out
fig = plt.figure(figsize=plt.figaspect(0.5))
ax = fig.add_subplot(1, 1, 1, projection='3d')
ax.set_title('3-D view of surface elevation')
# plot a 3D surface like in the example mplot3d/surface3d_demo
surf = ax.plot_surface(gwf.modelgrid.xcellcenters, 
                       gwf.modelgrid.ycellcenters, 
                       surface_elevation_cell_rl_ay_m,
                       rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
fig.colorbar(surf, shrink=0.5, aspect=10)
dis.top = surface_elevation_cell_rl_ay_m


#get_ipython().run_line_magic('matplotlib', 'inline')
#%% create observation package
obs = flopy.mf6.ModflowUtlobs(
    model  = gwf , # groundwater flow package to be used
    digits = 10, # default digits to print out
    print_input = True, 
    continuous= [['head_SA2', 'HEAD', point_centre_pond.lrc_loc],
                 ['head_SA3', 'HEAD', point_SA3.lrc_loc ],
                 ['head_SA4', 'HEAD', point_SA4.lrc_loc ]
                 ],
    pname = 'obs',
    filename='{}.obs'.format(fModName)
)


# find the cell which is intersected by the well

#wel_loc = (2, wel_id['cellids'][0][0], wel_id['cellids'][0][1])
#wel_spd = {1:[2, wel_id['cellids'][0][0], wel_id['cellids'][0][1], dch_wel, 0, 'dch_well']} # ordered dictionary to specify data for stress period 1

# wel = flopy.mf6.ModflowGwfwel(
#     gwf,
#     stress_period_data=wel_spd,
#     pname = 'wel'
# )

#np.ma.masked_equal(dis.sr.ygrid,1000) # very useful command to find specific file locations
#modelmap.sr.vertices
#flopy.plot.plotutil.cell_value_points
#
#
#modeLx_msect = flopy.plot.Modelc_mrossSection(model=mf, line={'Row': 0})



# oc = flopy.modflow.ModflowOc(mf, 
#                              stress_period_data=stress_period_data,
#                              compact=True)
# recharge package
#recharge_rate_mPday = {"0":0.1,"1":0.1}   #0.1 #0.01

# nrchop (int) – is the recharge option code. 
# 1: Recharge to top grid layer only 
# 2: Recharge to layer defined in irch 
# 3: Recharge to highest active cell (default is 3).
#nrchop = 1


# rch=flopy.modflow.ModflowRch(mf,
#                              rech   = recharge_rate_mPday,
#                              nrchop = nrchop,
#                              ipakcb = 1,
#                              stress_period_data   = {'0': adpoly_intersect.nodenumber , 
#                                        '1': adpoly_intersect.nodenumber}
#                              )
# rch = flopy.mf6.ModflowGwfrcha(
#     gwf, 
#     recharge   = 0.1,
#     #nrchop = 1
#     )


# rch_spd = [[None]]*(N**2-inactive_cells.size)
# count = 0
# for i in range(N):
#     for j in range(N):
#         if idomain[0, i, j] == idomain_rch:
#             rch_spd[count] = [0, i, j, 0.01, 0, 'rch_rain']
#             count += 1


#%% plot grid and idomain to show results plot the vertical view of the model
fig = plt.figure(figsize=(12, 9))
ax  = fig.add_subplot(1, 1, 1, aspect="equal")
modeLx_msect = flopy.plot.PlotCrossSection(model=gwf, line={'Row': 10})
arr = modeLx_msect.plot_array(dis.idomain[0,:,:])
modeLx_msect.plot_grid()
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)

#obs_ = mm.plot_bc("OBS")
ax.set_xlabel('X (m)', fontsize=40)
ax.set_ylabel('Y (m)', fontsize=40)
ax.set_title('idomain', fontsize=50)
#ax.colorbar(shrink=0.5, ax=ax)
#plt.colorbar(cax=ax)
cbar=plt.colorbar(arr, shrink=0.2, ax=ax)
cbar.ax.tick_params(labelsize=30)

#%%
fig = plt.figure(figsize=(15, 15))
ax = fig.add_subplot(1, 1, 1, aspect='equal')
#arr=g.plot(ax, a=dis.idomain, masked_values=[0], edgecolor='none', cmap='jet')
mm = flopy.plot.PlotMapView(model=gwf)
#mm.plot_ibound()
#flopy.plot.plot_shapefile(rf2shp, ax=ax, facecolor='yellow', alpha=0.25)

rch_ = mm.plot_bc("RCH")
chd_ = mm.plot_bc("CHD")
#obs_ = mm.plot_bc("OBS")

plt.plot(point_SA3.x,point_SA3.y,'ro',markersize=20)
plt.plot(point_centre_pond.x,point_centre_pond.y,'bo',markersize=20)
plt.plot(point_SA4.x,point_SA4.y,'go',markersize=20)

mm.plot_grid()
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
ax.set_xlabel('X (m)', fontsize=40)
ax.set_ylabel('Y (m)', fontsize=40)
ax.set_title('IDOMAIN', fontsize=50)
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
#import flopy.utils.binaryfile as bf
# Create the headfile and budget file objects
headobj       = flopy.utils.HeadFile(os.path.join(ws,fModName+'.hds'))
times_headobj = headobj.get_times()
cbcobj = flopy.utils.CellBudgetFile(os.path.join(ws,fModName+'.cbc'))
times_cbcobj = cbcobj.get_times()
#tsList = budObj.get_kstpkper() work but not needed.


# %%
# a list of times where the head distribution over time will be plotted.
times_output_list_day = [11, 101.0, 150, 200]
# find out the index of the head output (totim)
# that is cloest to the listed time for 
# plotting.  
index_of_totim_output_list_day = \
    [( abs(i-gwf.modeltime.totim) ).argmin() for i in times_output_list_day]

fig = plt.figure(figsize=(8, 3))
ax = fig.add_subplot(1, 1, 1)
#modeLx_msect = flopy.plot.Modelc_mrossSection(model=mf, line={'Column': 5})  
# this will onLy_m work when nrow is more than 1
##CM modeLx_msect = flopy.plot.Modelc_mrossSection(model=mf, line={'Row': 0})
modeLx_msect   = flopy.plot.PlotCrossSection(model=gwf, line={'Row': int(nrow/2) })
patches        = modeLx_msect.plot_ibound()
linecollection = modeLx_msect.plot_grid()
t = ax.set_title('Row 0 Cross-Section with idomain Boundary Conditions')
head = headobj.get_data(totim=times_output_list_day[1])
ax.plot(gwf.modelgrid.xycenters[0],head[-1,int(nrow/2),:])   
# head = headobj.get_data(totim=times_output_list_day[1])
# head[head==1e-30] = np.nan


# %% Plot the Head distribution at the end of stress period 1
fig = plt.figure(figsize=(8, 3))
ax  = fig.add_subplot(1, 1, 1)
# t   = ax.set_title('Head distribution at the end of stress period 1, day %i' %(times_output_list_day[0]))
# head = headobj.get_data(totim=times_output_list_day[0])
# modeLx_msect = flopy.plot.PlotCrossSection(model=gwf, line={'Row': 0})
# arr = modeLx_msect.plot_array(head)
# grd = modeLx_msect.plot_grid()
# ax.plot(gwf.modelgrid.xycenters[0] , head[-1,0,:], linewidth=5.0)
# plt.colorbar(arr, shrink=1, ax=ax)
times_head_ay = gwf.output.head().get_times()
head = gwf.output.head().get_data(totim = times_head_ay[-1])
times_bud_ay = gwf.output.head().get_times
bud = gwf.output.budget()

spdis = bud.get_data(text='DATA-SPDIS')[0]
qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)
pmv = flopy.plot.PlotMapView(gwf)
pmv.plot_array(head)
pmv.plot_grid(colors='white')
pmv.plot_vector(qx, qy, normalize=True, color="white")
pmv.contour_array(head, levels=[-5,-4,-3,-2], linewidths=3.)
ax.set_xlabel('X (m)')
ax.set_ylabel('Z (m)')


# %% Plot the Head distribution at the end of stress period 1
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
       

# %% obtain observation data from model output
# this may be replaced by pandas
time_array_obs_output_m, \
    hydrualic_head_SA2_time_array_m, \
    hydrualic_head_SA3_time_array_m, \
    hydrualic_head_SA4_time_array_m = \
    np.genfromtxt(ws+r'/0', 
    skip_header=1, 
    delimiter=',').T

#%% plot multiple graph to show changes of the results
fig, axes = plt.subplots(
    ncols=2,
    nrows=3,
    sharex=False,
    figsize=(9.3, 6.3),
    constrained_layout=True,
)

title_str ='kh = {:1.1e}'.format(hk_mPday) + '_ nlay = {:1.1e}'.format(nlay) \
    +'_ sy = {:1.1e}'.format(sy) \
#    + '_ nstp = {:1.1e}'.format(nstp) 
    
fig.suptitle(title_str,fontsize=10)

ax = axes[0,0]
ax.plot(field_data_df['time_elapsed_days'],
         field_data_df['sa2_watertable_rise_mm'] * mPmm,
         'bo',
         markevery=1000,
         label="SA2 Field")
ax.plot(time_array_obs_output_m,
        hydrualic_head_SA2_time_array_m - strt_m,
        'b-',
        label='SA2 modelled')
ax.set_xlabel('Time [s]')
ax.set_ylabel('Pressure head[m]')
ax.grid()
ax.set_ylim(-0.1,3)
ax.legend(loc="upper right",fontsize=10)

ax = axes[0,1]
ax.plot(field_data_df['time_elapsed_days'],
        field_data_df['sa3_watertable_rise_mm'] * mPmm,
        'ro',
        markevery=1000,
        label="SA3 Field")
ax.plot(time_array_obs_output_m,
        hydrualic_head_SA3_time_array_m - strt_m,
        'r-',
        label='SA3 modelled')
ax.grid()
ax.set_ylim(-0.1,3)
ax.set_xlabel('Time [s]')
ax.set_ylabel('Pressure head[m]')
ax.legend(loc="upper right",fontsize=10)

ax = axes[1,0]
ax.plot(time_array_obs_output_m,
        hydrualic_head_SA4_time_array_m - strt_m,
        'r-',
        label='SA4 modelled')
ax.plot(field_data_df['time_elapsed_days'],
        field_data_df['sa1_watertable_rise_mm'] * mPmm,
        'bo',
        markevery=1000,
        label="SA4 Field")
# ax.plot(time_array_obs_output_m,
#         hydrualic_head_SA2_time_array_m - strt_m,
#         'r.',
#         label='SA3 Measurement')
ax.grid()
ax.set_ylim(-0.1,3)
ax.set_xlabel('Time [s]')
ax.set_ylabel('Pressure head[m]')
ax.legend(loc="upper right",fontsize=10)

ax = axes[1,1]
ax.plot(field_data_df['time_elapsed_days'],
        field_data_df['surface_water_depth_m'],
        'ko',
        markevery=1000,
        label="measured")
# ax.plot(time_array_obs_output_m,
#         hydrualic_head_SA2_time_array_m - strt_m,
#         'r.',
#         label='SA3 Measurement')
ax.grid()
ax.set_ylim(-0.1,1)
ax.set_xlabel('Time [s]')
ax.set_ylabel('Water depth [m]')
ax.legend(loc="upper right",fontsize=10)


ax =  axes[2,1]
modeLx_msect   = flopy.plot.PlotCrossSection(model=gwf, 
                                             line={'Row': int(nrow/2) })
patches        = modeLx_msect.plot_ibound()
linecollection = modeLx_msect.plot_grid()
t = ax.set_title('Head a Cross-Section')
#head = headobj.get_data(totim=times_output_list_day[1])
for i in index_of_totim_output_list_day :
    head = gwf.output.head().get_data(totim = times_head_ay[i])
    ax.plot(gwf.modelgrid.xycenters[0],
            head[-1,int(nrow/2),:],
            label ='day ' +str(times_head_ay[i]) 
            )   
ax.legend(loc="lower right",fontsize=10)



ax = axes[2,0]
ax.plot(stress_period_end_time_days_ay,
        recharge_rate_spd_ay,
        'ko',
        label="measuched")
# ax.plot(time_array_obs_output_m,
#         hydrualic_head_SA2_time_array_m - strt_m,
#         'r.',
#         label='SA3 Measurement')
ax.grid()
#ax.set_ylim(-0.1,1)
ax.set_xlabel('Time [s]')
ax.set_ylabel('Recharge [m]')
ax.legend(loc="upper right",fontsize=10)



plt.show()


fname_save=title_str.replace(';', ' ').replace('+', '').replace('e-', 'ne') \
    .replace('=', '_').replace(' ', '')
print(fname_save)
plt.savefig(fname_save+'.png',dpi=300)


