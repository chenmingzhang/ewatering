# -*- coding: utf-8 -*-
"""

@author:  Chenming Zhang, Yuanchi Xiao

procedure to run the code:
    
    1. In the lab PC environment, open up Laboratory Software -> Anaconda3 (64-bit) -> Anaconda Powershell Prompt (Anaconda3)
    2. Packages to be installed (may need to answer y during the installation): 
        conda install -c main rasterio=1.3.8
        conda install -c conda-forge gdal=3.6.2
        conda install -c conda-forge flopy
        conda install matplotlib=3.7.1
        #conda install osgeo
        pip install gdal ==3.6.2
        #conda update Pillow
       Avoid mixing use of flopy and conda interchangeablly that may cause dependency issues
       
       for python3.12.x
       conda install -c conda-forge gdal
       conda install matplotlib
       
    3. enter spyder in the Anaconda Powershell Prompt
        
as of 241001
3.12.3 | packaged by conda-forge | (main, Apr 15 2024, 18:20:11) [MSC v.1938 64 bit (AMD64)]
numpy version: 1.26.4
matplotlib version: 3.9.2
flopy version: 3.8.1


#CZ241001 matplotlib qt needs to be executed to ensure the contour map is shown 
on top of a map.

%matplotlib qt
or 
import matplotlib as mpl
mpl.use('TkAgg')


Stress period 1 

Stress period 1 (30 years): boundary + river, to simulate the conditions before farming
Stress period 2 (2 years): boundary + river + farming, for 50 years, to simulate the additional flow and salt load to the river caused by irrigation
Stress period 3 (3 years): boundary + river + farming + pump, to simulate how effective salt interception scheme can reduce salt load
Stress period 4: boundary + river + farming + pump + watering, to simulate how much additional salt is introduced to the river due to watering


"""
# %%## 
# =============================================================================
# ===== CIVL4145 Course Project 
# ===== USING MODFLOW 6 TO SIMULATE GROUNDWATER FLOW AT A FLOODPLAIN OF The RIVER MURRAY
# =============================================================================

# =============================================================================
# ===== 0.0 PYTHON ENVIRONMENT SET-UP =========================================
# =============================================================================
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import flopy
import numpy as np
from osgeo import gdal
import pandas as pd
import matplotlib as mpl
import sys
import os

print(sys.version)
print(f"numpy version: {np.__version__}")
print(f"matplotlib version: {mpl.__version__}")
print(f"flopy version: {flopy.__version__}")

# Specifying executable location, working folder and simulation name
#mf6exe = '//puffball.labs.eait.uq.edu.au/uqczhan2/Desktop/modflow6/ewatering_flopy/mf6.exe'
mf6exe =  os.path.join(os.getcwd(),'mf6.exe') #'D:/OneDrive/OneDrive - The University of Queensland/Teaching/CIVL4145/2024/ewatering_flopy/mf6.exe'
ws = '' # working directory
sim_name      = "5_stress_pd"
layer_1_top   = "Coonambidgal_dem.tif"
layer_2_top   = "Monoman_dem.tif"
recharge_data_irrigation_only = pd.read_csv("recharge_irrigation_only.csv")
recharge_data = pd.read_csv("recharge_irrigation_inundation.csv")
river_data    = pd.read_csv("river.csv")
sis_well_data = pd.read_csv("sis_well.csv")

layer_1_top_data                = gdal.Open(layer_1_top)
#layer_1_top_band                = layer_1_top_data.GetRasterBand(1)
#layer_1_top_nodataval           = layer_1_top_band.GetNoDataValue()
layer_1_top_m_mtx               = layer_1_top_data.ReadAsArray().astype(np.float64)
#top = layer_1_top_m_mtx
layer_2_top_data                = gdal.Open(layer_2_top)
#layer_2_top_band                = layer_2_top_data.GetRasterBand(1)
#layer_2_top_nodataval           = layer_2_top_band.GetNoDataValue()
layer_2_top_m_mtx               = layer_2_top_data.ReadAsArray().astype(np.float64)
# a mask matrix to identify the highland where the elevation is greater than 20 mAHD
mask_highland                   = layer_1_top_m_mtx >= 20


#defining values for plotting
fs_title_1 = 12
fs_title_sub = 10
mpl.rcParams['image.cmap'] = 'jet'
# %%# --- 
# =============================================================================
# ===== 1.0 DEFINING BASIC FloPy MODEL PROPERTIES =============================
# =============================================================================

# ===== 1.1 Setting length and time units
length_units = "meters"
time_units = "days"
no_days_in_a_year = 365
# Setting three stress periods perlen(10000 days, 351 days, 100 days)
#perlen = [1*10000, no_days_in_a_year*2,no_days_in_a_year*10,no_days_in_a_year*1,no_days_in_a_year*2]  # Simulation time [days]
perlen = no_days_in_a_year* np.array([30, 2,3,1,2])
perlen_cum = np.cumsum(perlen)
nper = len(perlen)
#nstp = [1, 10, no_days_in_a_year*3, no_days_in_a_year,no_days_in_a_year/5]  # orginal
#nstp = [1, 10, int(no_days_in_a_year/10), int(no_days_in_a_year/10), int(no_days_in_a_year/10)]  # working 

nstp = [1, 10, int(no_days_in_a_year/20), int(no_days_in_a_year/20), int(no_days_in_a_year/20)]   

tsmult = [1.0, 1.0, 1.0, 1.0,1.0]

# ===== 1.2 Setting model geometry
Lx = 7400    # Length of domain in x-direction
Ly = 5800    # Length of domain in y-direction
delr = 40.0  # Row width (height) [m] (i.e., dy)
delc = 40.0  # Column width [m] (i.e., dx)
# grid location info
xoff = 482460 # lower left x-coordinate in UTM meters
yoff = 6230437 # lower left y-coordinate
nrow = int(np.round(Ly/delr)) # Number of rows
# Add extra 2 columns for west & east boundary cells (affects x plot locations)
ncol = int(np.round(Lx/delc)) # Number of columns
nlay = 2  # Number of layers, only coonambidgal and Monoman sands are included
#delz = [3.0, 7.0, 2.0]  # Layer thickness [m] (starting at top layer),need to use dem 
Coonambidgal_thickness = layer_1_top_m_mtx - layer_2_top_m_mtx
#top = 12.0  # Elevation of top of the model [m], need to chek dem
basez = -26.0  # Elevation of base of the model [m]
# Settng elevation of bottom of each layer
# botz = []
# botz.append(top-delz[0])
# for k in range(1,nlay):
#     botz.append(botz[k-1]-delz[k])

botm = np.zeros((nlay, nrow, ncol), dtype=float)
# for k in np.arange(nlay):
botm[0,:,:] = layer_2_top_m_mtx #coonambidgal clay bottom, also the top of the monoman sand 
botm[1,:,:] = basez*np.ones((nrow, ncol), dtype=float) #sand bottom, also the top of the monoman sand 

# %%===== 1.3 Setting hydraulic and solute transport properties
# Porosity [-] Ss and Sy, essential:note that porosity is only effctive for gwt
#prsity = 0.30  # Porosity
nclay = 0.03    # Porosity of clay
nsand = 0.15    # Porosity of sand
prsity = np.ones((nlay, nrow, ncol), dtype=float) # Establish variable
prsity[0,:,:] = nclay   # Set layer 2 porosity to nsand
prsity[1,:,:] = nsand   # Set layer 0 porosity to nsand  
ss_sand = 1e-4
ss_clay = 1e-4
ss = np.ones((nlay, nrow, ncol), dtype=float) # Establish variable
ss[0,:,:] = ss_clay
ss[1,:,:] = ss_sand

sy_sand = 0.15  # specific yield of sand
sy_clay = 0.03  # specific yield of clay
sy = np.ones((nlay, nrow, ncol), dtype=float)
sy[0,:,:] = sy_clay
sy[1,:,:] = sy_sand

# Hydraulic conductivity [m/d]
# k11 = 44.0   # Horizontal hydraulic conductivity [m/d]
kclay = 0.05 # Horizontal hydraulic conductivity of clay [m/d]
# ksand_floodplain = 1.0 # Horizontal hydraulic conductivity of sand [m/d]
# ksand_highland   = 8.0 # loxton sand
ksand_floodplain = 0.5 # Horizontal hydraulic conductivity of sand [m/d]
ksand_highland   = 4 # loxton sand
k11 = np.ones((nlay, nrow, ncol), dtype=float) # Establish variable
k11[0,:,:] = kclay
k11_2      = ksand_floodplain*np.ones((nrow, ncol))
k11_2[mask_highland] = ksand_highland
k11[1,:,:] = k11_2

k33 = 0.01*k11.copy()  # Vertical hydraulic conductivity [m/d]


# %% plot the hydraulic conductivity
fig = plt.figure(dpi=200)
ax1 = plt.subplot(1,2,1)
ax1.set_aspect('equal')
im1=ax1.imshow(k11[0,:,:]) 
ax1.set_xlabel('Distance (m)')
ax1.set_ylabel('Distance (m)')
ax1.set_title('Layer 1 ',fontsize=fs_title_sub)
#im1.set_clim(0,5)
ax2 = plt.subplot(1,2,2)
im2=ax2.imshow(k11[1,:,:]) 
ax2.set_aspect('equal')
ax2.set_xlabel('Distance (m)')
ax2.set_ylabel('Distance (m)')
ax2.set_title('Layer 2',fontsize=fs_title_sub)
#im2.set_clim(0,5)
#plt.tight_layout()
#cbar_ax = fig.add_axes([0.35, 0.18, 0.4, 0.02])
# cbar = fig.colorbar(im1, cax=cbar_ax,location='bottom',orientation='horizontal')
# cbar.ax.tick_params(axis='x', direction='in')
plt.figtext(0.3, 0.8, 'Hydraulic conductiivity (m/day)',fontsize=fs_title_1)
plt.savefig('Hydraulic_conductivity.png')

# %% Dispersivity (nitrate) inactivated
# al = 10              # Longitudinal dispersivity [m]
# trpt = 0.1             # Ratio of transverse to longitudinal dispersitivity
# trpv = 0.1           # Ratio of vertical to longitudinal dispersitivity
# ath1 = al * trpt       # Transverse dispersivity [m]
# ath2 = al * trpv       # Vertical dispersivity [m]
# sconc = 35.0            # Starting soute concentration (C(t=0))
# dmcoef = 1.468e-4     # Molecular diffusion coefficient [m**2/d]

# Recharge
# rech = 0.012 # Recharge rate [m/d], 12 mm/day

# %% Specify saturated thickness method
# (0=constant (confined?), >0=varies with head in cell (unconfined), <0=variable)
icelltype = 1

# Advection solution options
mixelm = -1
# mixelem = 0, the standard finite-difference method with upstream
#              or central-in-space weighting, depending on the value of NADVFD;
# mixelem = 1, the forward-tracking method of characteristics (MOC);
# mixelem = 2, the backward-tracking modified method of characteristics (MMOC);
# mixelem = 3, the hybrid method of characteristics (HMOC) with MOC or MMOC
#              automatically and dynamically selected;
# mixelem = -1, the third-order TVD scheme (ULTIMATE).

# ===== 1.4 Set initial conditions
h = 16.0 # Hydraulic head on western boundary
# h2 = 9.6 # Hydraulic head on eastern boundary
# havg = (h1+h2)/2 # Average hydraulic head to use as initial cond
strt = np.ones((nlay, nrow, ncol), dtype=float)*h
#    l, r, c 

# Exlude cells from calculation using idomain flag
# idomain = 0 --> cells does not exist in simulation, but values are written to it
# idomain = 1 --> cell exists in simulation
# idomain = -1 --> cell does not exist in simulation, but flow passes through it
idomain = np.ones((nlay, nrow, ncol), dtype=int) 


well   = []
welspd = {}
for i in range(len(sis_well_data)):
#    well.append([(sis_well_data.l[i]-1,sis_well_data.r[i]-1,sis_well_data.c[i]-1),\
#            -sis_well_data.Q[i],sis_well_data.concentration[i]])
    well.append([ (sis_well_data.l[i]-1,sis_well_data.r[i]-1,sis_well_data.c[i]-1),\
            -sis_well_data.Q[i] ,sis_well_data.concentration[i]])
        
welspd = {0:None,1:None,2:well,3:well,4:well}    # CZ241001 they do not like None


# flopy.mf6.ModflowGwfwel(gwf,
#                         print_input=True,
#                         print_flows=True,
#                         stress_period_data=welspd,
#                         save_flows=True,
#                         auxiliary="WELL",
#                         pname="WEL-1",
#                         filename="{}.wel".format(gwfname),)
# %%===== 1.5 Define recharge information
# Recharge stress perdiod data needed if recharge method 1 (general recharge
# package is used).
rchspd = []
rch_inundation_irrigation = []
rch_inundation_only       = []
for i in range(len(recharge_data)):
    row    = recharge_data.row[i]-1
    column = recharge_data.column[i]-1
    rech   = recharge_data.Rech[i]
    conc   = recharge_data.Conc[i]
        #             [(lay, row, col), recharge, conc]
    rch_inundation_irrigation.append([(  0,   row,   column),     rech,  0.001])

for i in range(len(recharge_data_irrigation_only)):
    row    = recharge_data_irrigation_only.row[i]-1
    column = recharge_data_irrigation_only.column[i]-1
    rech   = recharge_data_irrigation_only.Rech[i]
    conc   = recharge_data_irrigation_only.Conc[i]
        #             [(lay, row, col), recharge, conc]
    rch_inundation_only.append([(  0,   row,   column),     rech,  0.001])
    
rchspd = {0: None,1: rch_inundation_only,2:rch_inundation_only,3:rch_inundation_irrigation, 4: rch_inundation_only}

# ===== 1.6 Define general head boundaries
left_ghb_m  = 17.6
right_ghb_m = 18.0

ghbspd = []
ghb = []
for col in range(ncol):
    stage = left_ghb_m + (col / (ncol - 1)) * (right_ghb_m - left_ghb_m)  # Linearly increasing stage
    #ghb.append([1, nrow-1, col, stage, 1000])  # [layer, row, column, stage, conductance]
    ghb.append([ (1,nrow-1, col),stage, 1000,layer_1_top_m_mtx[nrow-1,col] ]) #,  1000,stage])  # [(k, 0, ncol - 1), top, ghbcond, 35.0] 
    # https://github.com/MODFLOW-USGS/modflow6-examples/blob/c5b10d3bcd87c06810dc954a788f7fb739ca2cc2/scripts/ex-gwt-henry.py#L137
ghbspd = {0: ghb, 1: ghb,2: ghb,3:ghb,4:ghb}


# %% ===== 1.7 Define river boundaries
 
riv    = []
rivspd = []
for i in range(len(river_data)):
    row    = river_data.row[i]-1
    column = river_data.column[i]-1
    stage  = river_data.Stage[i]
    conc   = river_data.Conc[i]
    conductance = river_data.Conductanc[i]
    if river_data.Bottom[i] < layer_1_top_m_mtx[row,column]:
        botm[0,row,column]  = river_data.Bottom[i] - 0.2
        rbot = river_data.Bottom[i]
    else:
        rbot = river_data.Bottom[i]
    # CZ241001 this seems to be revised in flopy
    riv.append([(0, row, column), stage, conductance, rbot,0])
rivspd = {0: riv, 1: riv, 2:riv,3:riv,4:riv}
        
# ===== 1.7 Define solver settings
nouter, ninner = 100, 300
hclose, rclose, relax = 1e-6, 1e-6, 1.0
percel = 1.0  # HMOC parameters
itrack = 3
wd = 0.5
dceps = 1.0e-5
nplane = 1
npl = 0
nph = 16
npmin = 2
npmax = 32
dchmoc = 1.0e-3
nlsink = nplane
npsink = nph

# ====== Set static temporal data used by tdis file
tdis_rc = []
tdis_rc.append((perlen, nstp, 1.0))

# %%=============================================================================
# ===== 2.0 CREATE FLOW MODEL OBJECTS AND DEFINE FLOW PACKAGES ================
# =============================================================================
name = "Ewatering"
gwfname = "gwf-" + name
sim_ws = ws
sim = flopy.mf6.MFSimulation(sim_name=sim_name,
                             sim_ws=sim_ws,
                             exe_name=mf6exe)

# ===== 2.1 Defining MODFLOW 6 time discretization
tdis_rc = []
for i in range(nper):
    tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

flopy.mf6.ModflowTdis(sim,
                      nper=nper,
                      perioddata=tdis_rc,
                      time_units=time_units)
        
# ===== 2.2 Defining MODFLOW 6 groundwater flow model
gwf = flopy.mf6.ModflowGwf(sim,
                           modelname=gwfname,
                           save_flows=True,
                           model_nam_file="{}.nam".format(gwfname),)             
# ===== 2.3 Defining MODFLOW 6 solver for flow model
imsgwf = flopy.mf6.ModflowIms(sim,
                              print_option="SUMMARY",
                              outer_dvclose=hclose,
                              outer_maximum=nouter,
                              under_relaxation="NONE",
                              inner_maximum=ninner,
                              inner_dvclose=hclose,
                              rcloserecord=rclose,
                              linear_acceleration="BICGSTAB",
                              scaling_method="NONE",
                              reordering_method="NONE",
                              relaxation_factor=relax,
                              filename="{}.ims".format(gwfname),)

sim.register_ims_package(imsgwf, [gwf.name])

# ===== 2.4 Defining MODFLOW 6 discretization package
flopy.mf6.ModflowGwfdis(gwf,
                        length_units=length_units,
                        nlay=nlay,
                        nrow=nrow,
                        ncol=ncol,
                        delr=delr,
                        delc=delc,
                        xorigin=xoff,
                        yorigin=yoff,
                        top=layer_1_top_m_mtx,
                        botm=botm,
                        idomain=idomain,
                        filename="{}.dis".format(gwfname),)

# ===== 2.5 Defining MODFLOW 6 node-property flow package
flopy.mf6.ModflowGwfnpf(gwf,
                        save_flows=False,
                        icelltype=icelltype,
                        k=k11,
                        k33=k33,
                        save_specific_discharge=True,
                        filename="{}.npf".format(gwfname),)

# ===== 2.6 Defining MODFLOW 6 initial conditions package for flow model
flopy.mf6.ModflowGwfic(gwf,
                       strt=strt,
                       filename="{}.ic".format(gwfname),)


# ===== 2.7 Define MODFLOW 6 storage package
sto = flopy.mf6.ModflowGwfsto(gwf,
                              ss=ss,
                              sy=sy,
                              filename="{}.sto".format(gwfname),)


#ghbspd = {0: ghb, 1: ghb,2: ghb,3:ghb ,4:ghb}
# ===== 2.8 Defining MODFLOW 6 constant head package
flopy.mf6.ModflowGwfghb(gwf,
                        maxbound=len(ghbspd),
                        stress_period_data=ghbspd,
                        save_flows=True,
                        auxiliary="GHB",
                        pname="GHB-1",
                        filename="{}.ghb".format(gwfname),)

# ===== 2.9 Defining MODFLOW well package
flopy.mf6.ModflowGwfwel(gwf,
                        print_input=True,
                        print_flows=True,
                        stress_period_data=welspd,
                        save_flows=True,
                        auxiliary="WELL",
                        pname="WEL-1",
                        filename="{}.wel".format(gwfname),)

# ===== 2.10 Defining MODFLOW recharge package
# Method 1 - Using general recharge package
flopy.mf6.ModflowGwfrch(gwf,
                        print_input=True,
                        print_flows=True,
                        stress_period_data=rchspd,
                        save_flows=True,
                        auxiliary="RCH",
                        pname="RCH-1",
                        filename="{}.rch".format(gwfname),)

flopy.mf6.ModflowGwfriv(gwf,
                        print_input=True,
                        print_flows=True,
                        stress_period_data=rivspd,
                        save_flows=True,
                        auxiliary="RIV",
                        pname="RIV-1",
                        filename="{}.riv".format(gwfname),)

# Method 2 - Using array-based recharge package
# flopy.mf6.ModflowGwfrcha(gwf,
#                         print_input=True,
#                         print_flows=True,
#                         recharge = rech,
#                         save_flows=False,
#                         auxiliary="CONCENTRATION",
#                         pname="RCH-1",
#                         filename="{}.rch".format(gwfname),)

# ===== 2.11Defining MODFLOW 6 output control package for flow model
flopy.mf6.ModflowGwfoc(gwf,
                       head_filerecord="{}.hds".format(gwfname),
                       budget_filerecord="{}.cbb".format(gwfname),
                       
                       headprintrecord=[("COLUMNS", 10, "WIDTH", 15,
                                         "DIGITS", 6, "GENERAL")],
                       budgetcsv_filerecord="budget.csv",
                       saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
                       printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],)

# %%#=============================================================================
# #===== 3.0 CREATE TRANSPORT MODEL OBJECTS AND DEFINE TRANSPORT PACKAGES ======
# #=============================================================================
# gwtname = "gwt_" + name
# gwt = flopy.mf6.MFModel(sim,
#                         model_type="gwt6",
#                         modelname=gwtname,
#                         model_nam_file="{}.nam".format(gwtname),)
# gwt.name_file.save_flows = True

# # ===== 3.1 Create iterative model solution and register the gwt model with it
# imsgwt = flopy.mf6.ModflowIms(sim,
#                               print_option="SUMMARY",
#                               outer_dvclose=hclose,
#                               outer_maximum=nouter,
#                               under_relaxation="NONE",
#                               inner_maximum=ninner,
#                               inner_dvclose=hclose,
#                               rcloserecord=rclose,
#                               linear_acceleration="BICGSTAB",
#                               scaling_method="NONE",
#                               reordering_method="NONE",
#                               relaxation_factor=relax,
#                               filename="{}.ims".format(gwtname),)
# sim.register_ims_package(imsgwt, [gwt.name])

# # ===== 3.2 Defining MODFLOW 6 transport discretization package
# flopy.mf6.ModflowGwtdis(gwt,
#                         nlay=nlay,
#                         nrow=nrow,
#                         ncol=ncol,
#                         delr=delr,
#                         delc=delc,
#                         top=top,
#                         botm=botm,
#                         idomain=1,
#                         filename="{}.dis".format(gwtname),)

# # ===== 3.3 Defining MODFLOW 6 transport initial concentrations

# iconc = np.zeros((nlay, nrow, ncol), dtype=float) #initial concentration
# # for k in np.arange(nlay):
# iconc[0,:,:] = 0*np.ones((nrow, ncol), dtype=float) #coonambidgal clay salt concentration
# iconc[1,:,:] = sconc*np.ones((nrow, ncol), dtype=float) #monoman sand salt concentration

# flopy.mf6.ModflowGwtic(gwt,
#                         strt=iconc,
#                         filename="{}.ic".format(gwtname),)

# # ===== 3.4 Defining MODFLOW 6 transport advection package
# if mixelm == 0:
#     scheme = "UPSTREAM"
# elif mixelm == -1:
#     scheme = "TVD"
# else:
#     raise Exception()

# flopy.mf6.ModflowGwtadv(gwt,
#                         scheme=scheme,
#                         filename="{}.adv".format(gwtname),)


# # ===== 3.5 Defining MODFLOW 6 transport dispersion package
# if al != 0:
#     flopy.mf6.ModflowGwtdsp(gwt,
#                             xt3d_off=True,
#                             alh=al,
#                             ath1=ath1,
#                             ath2=ath2,
#                             diffc=dmcoef,
#                             filename="{}.dsp".format(gwtname),)

# #===== 3.6 Defining MODFLOW 6 transport mass storage package
# #      (formerly "reaction" package in MT3DMS)
# flopy.mf6.ModflowGwtmst(gwt,
#                         porosity=prsity,
#                         first_order_decay=False,
#                         decay=None,
#                         decay_sorbed=None,
#                         sorption=None,
#                         bulk_density=None,
#                         distcoef=None,
#                         filename="{}.mst".format(gwtname),)

# ===== 3.7 Defining MODFLOW 6 transport constant concentration package
# flopy.mf6.ModflowGwtcnc(gwt,
#                         maxbound=len(cncspd),
#                         stress_period_data=cncspd,
#                         save_flows=False,
#                         pname="CNC-1",
#                         filename="{}.cnc".format(gwtname),)

# # ===== 3.8 Defining MODFLOW 6 transport source-sink mixing package
# sourcerecarray = [("GHB-1", "AUX", "GHB"),
#                   ("RCH-1", "AUX", "RCH")]
# flopy.mf6.ModflowGwtssm(gwt,
#                         sources=sourcerecarray,
#                         filename="{}.ssm".format(gwtname),)

# # ===== 3.9 Defining MODFLOW 6 transport output control package
# flopy.mf6.ModflowGwtoc(gwt,
#                         budget_filerecord="{}.cbc".format(gwtname),
#                         concentration_filerecord="{}.ucn".format(gwtname),
#                         concentrationprintrecord=[
#                             ("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6,
#                             "GENERAL")],
#                         saverecord=[("CONCENTRATION", "ALL"),
#                                     ("BUDGET", "LAST")],
#                         printrecord=[("CONCENTRATION", "LAST"),
#                                     ("BUDGET", "LAST")],)

# # Defining MODFLOW 6 flow-transport exchange mechanism
# flopy.mf6.ModflowGwfgwt(sim,
#                         exgtype="GWF6-GWT6",
#                         exgmnamea=gwfname,
#                         exgmnameb=gwtname,
#                         filename="{}.gwfgwt".format(name),)

#%%=============================================================================
# ===== 4.0 CREATE MODFLOW 6 INPUT FILES AND RUN THE SIMULATION ===============
# =============================================================================

# ===== 4.1 Write input files
sim.write_simulation()

# ===== 4.2 Run the simulation
success, buff = sim.run_simulation()
assert success, "MODFLOW 6 did not terminate normally."

# %%    ----

# %%=============================================================================
# ===== 5.0 POST-PROCESS SIMULATION RESULTS ===================================
# =============================================================================
import matplotlib.image as image
import rasterio
from rasterio.plot import show
# ===== 5.1 Extracting simulation data ========================================
# Get the MF6 flow model output
# ===== Output method
headobj_mf6_head = gwf.output.head()
head = headobj_mf6_head.get_alldata()
timesteps = headobj_mf6_head.times
base_map = rasterio.open('map_murtho_large.tif')
map_band1 = base_map.read(1)
height, width = map_band1.shape
cols, rows = np.meshgrid(np.arange(width), np.arange(height)) 
xs, ys = rasterio.transform.xy(base_map.transform, rows, cols) 
xcoords = np.array(xs)
ycoords = np.array(ys)
#base_map = base_map.band(1).shape[0]
# timesteps_h = ucnobj_mf6_head.times
# head_end = head[-1]
# ===== Binary file method
# fname = 'Run000/gwf-Base.hds'
# hdobj = flopy.utils.HeadFile(fname)
# head = hdobj.get_data()

# ===== Get the MF6 transport model output (concentration)
# ucnobj_mf6_conc = gwt.output.concentration()
# conc = ucnobj_mf6_conc.get_alldata()
# timesteps = ucnobj_mf6_conc.times

# ===== Set timestep to use
targetTime_sp1 = np.array([0])
targetTime_sp2 = perlen[0]+np.array([200,300,400,500,600])
targetTime_sp3 = perlen[0]+perlen[1]+np.int32(no_days_in_a_year * np.linspace(0,3,10))
targetTime_sp4 = perlen[0]+perlen[1]+perlen[2]+np.int32(np.linspace(0,no_days_in_a_year,10))
targetTime_sp5 = perlen[0]+perlen[1]+perlen[2]+perlen[3]+np.int32(np.linspace(0,no_days_in_a_year*2,10))


targetTime_sp1 = np.array([0])
targetTime_sp2 = perlen[0]+np.int32( np.linspace(0,  perlen[1], 5 )  )
targetTime_sp3 = perlen[0]+perlen[1]+np.int32( np.linspace(0,  perlen[2], 5 )  )
targetTime_sp4 = perlen[0]+perlen[1]+perlen[2] +np.int32( np.linspace(0,  perlen[3], 5 )  )
targetTime_sp5 = perlen[0]+perlen[1]+perlen[2]++np.int32( np.linspace(0,  perlen[4], 5 )  )

targetTime = np.concatenate([targetTime_sp1,targetTime_sp2,targetTime_sp3,targetTime_sp4,targetTime_sp5])
timeIDX=np.zeros(np.size(targetTime),dtype='int32')
klay=1
for k in range(len(targetTime)):
    print(k)
    timeIDX[k] = np.int32([i for i, v in enumerate(timesteps) if v >= targetTime[k]][0])
    dispTime = timesteps[timeIDX[k]]
    # ===== Set layer to plot
    fig = plt.figure(figsize=(12, 10),dpi=150, tight_layout=True)
    ax1 = fig.add_subplot(1, 1, 1, aspect="equal")
    ax1.set_title("Hydrualic Head: Layer " + "{}".format(klay) +
             " (Time = " + "{:.1f}".format(dispTime/no_days_in_a_year) + " "
             "{}".format("years") + ")",
             loc='left')
    plt.xlabel("Distance along x-axis [m]")
    plt.ylabel("Distance along y-axis [m]")
    mapview = flopy.plot.PlotMapView(model=gwf,\
                                 extent=(np.min(xcoords),np.max(xcoords),\
                                         np.min(ycoords),np.max(ycoords)),layer=1)
    show(base_map)
    theLevels = np.linspace(16,18,num=6)
    contour_set = mapview.contour_array(head[timeIDX[k],klay,:,:],
                                    levels=theLevels,
                                    colors='k',
                                    linestyles='--')
    plt.clabel(contour_set, inline=1, fontsize=10)
    # plt.colorbar(contour_set, shrink=0.250)
    quadmesh = mapview.plot_array(head[timeIDX[k],klay,:,:], alpha=0.2,clim=(16,18))
    
    #cb = plt.colorbar(quadmesh, shrink=0.25)
    #linecollection=mapview.plot_grid()
    quadmesh = mapview.plot_bc("GHB",color='black')
    quadmesh = mapview.plot_bc("WEL-1",kper=3,color='red')
    mapview2 = flopy.plot.PlotMapView(model=gwf,\
                                 extent=(np.min(xcoords),np.max(xcoords),\
                                         np.min(ycoords),np.max(ycoords)),layer=0)
    quadmesh = mapview2.plot_bc("RIV",color='green')
    quadmesh = mapview2.plot_bc("RCH",kper=3,color='blue',alpha=0.2)
    # quadmesh3 = mapview.plot_bc("",color='red')
    # ax2 = fig.add_subplot(2, 1, 2, aspect="auto")
    # plt.bar(0,2)
    plt.savefig(f'{targetTime[k]}.png')
    plt.close()

# %%===== 5.4 Plotting cross-section of head and solute data
rowID = 70 # Specify the row to use for cross section
colID = 111
klay  = 1
fig = plt.figure(figsize=(8, 8),dpi=150, tight_layout=True)
# Hydraluic heads vertical section
ax = fig.add_subplot(1, 1, 1,aspect=30)
ax.set_title("Hydraulic head: Vertical section - Row " +
             "{}".format(rowID) + " (Time = " + "{:.1f}".format(dispTime/no_days_in_a_year  ) + " "
             "{}".format('years') + ")", loc='left')
plt.xlabel("Distance along x-axis [m]")
plt.ylabel("Elevation (z-axis) [m]")
xsect = flopy.plot.PlotCrossSection(model=gwf, line={"column": rowID})
csa = xsect.plot_array(head, head=head,alpha=0.5)
cb = plt.colorbar(csa, shrink=0.25)
contour_set = xsect.contour_array(head,
                                  head=head,
                                  levels=theLevels,
                                  colors="k",
                                  linestyles='--')
plt.clabel(contour_set, fmt="%.1f", colors="k", fontsize=11)
# patches = ysect.plot_bc("GHB",color='pink')
# patches = xsect.plot_bc("WEL-1",color='grey')
linecollection = xsect.plot_grid()
plt.tight_layout()


plt.figure()
plt.plot(np.array(timesteps)/no_days_in_a_year,head[:,klay,rowID,colID])
plt.xlabel('Time (year)')
plt.ylabel('Head (mAHD)')
plt.title('Head at SA2')
plt.grid()
plt.savefig('head_sa2.jpg')
plt.close('all')
# # create modpath files
# mpnamf = f"{sim_name}_mp_forward"

# # create basic forward tracking modpath simulation
# mp = flopy.modpath.Modpath7('mgwf-mp', flowmodel=gwf)
# mpsim = flopy.modpath.Modpath7Sim(mp,referencetime=[4,0,1])
# # write modpath datasets
# mpsim.write_file()
