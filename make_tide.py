'''
this is an example to generate a ROMS tide file for TPXO8-atlas tidal product. Till this moment I couldn't found
a good python script to do so, so this script is created for this purpose. TPXO8-atlas is a high resolution global
model of ocean tides. it is available here:

http://volkov.oce.orst.edu/tides/tpxo8_atlas.html

download the 8 major constitutes with 1/30 deg resotuion. Other constitutes are also supported, but haven't been
tested. Only netCDF format is supported in this version, so make sure you download the right files.

The script 'make_remap_weights_file.py' is used to generate weight file between two grids, in this example the
Glacier Bay and TPXO8 grid files. TPXO8 grid file is kept in the same directory as data files:

/Volumes/R1/scratch/chuning/data/tpxo8nc/grid_tpxo8atlas_30_v1.nc

Run this script before 'make_tide.py' to generate the weight files (remap_*.nc).

To use the remap function a new set of grid is created (CGrid_TPXO8). I followed the structure of GLORY CGrid 
from the PYROMS packaage itself, so that it is consistent with PYROMS.

---------------------------------------------------------------------------------------------------------------------
tidal parameters are converted to elliptical parameters using ap2ep, origionally written in MATLAB and later
translated to python package PyFVCOM. See

CGrid_TPXO8/tidal_ellipse.py for more details.

-Chuning Wang, Dec 18 2016

'''

import netCDF4 as nc
import numpy as np

import pyroms

grd1 = 'fjord_test'

# read ROMS grid
dst_grd = pyroms.grid.get_ROMS_grid(grd1)
x = dst_grd.hgrid.x_rho
y = dst_grd.hgrid.y_rho
msk = dst_grd.hgrid.mask_rho
h = dst_grd.vgrid.h
eta, xi = msk.shape

# -------------------------------------------------------------------------
# tidal constituents
# consts1 = ['Q1', 'O1', 'P1', 'K1', 'N2', 'M2', 'S2', 'K2']
# # consts2 = ['MF']
# consts2 = []
# consts = consts1+consts2
# # consts = ['M2'] 
# consts_num = len(consts)
# 
# # define tidal constituents names and periods
# tide_name = np.array([ list('Q1  '), list('O1  '), list('P1  '), list('K1  '),
#                        list('N2  '), list('M2  '), list('S2  '), list('K2  ')])
# #                        list('MF  ')])
# tide_period = np.array([26.8683567047119, 25.8193397521973, 24.0658893585205, 23.9344692230225,
#                         12.6583499908447, 12.420599937439, 12, 11.9672346115112])
# #                         13.66079*24])

consts = ['M2']
consts_num = len(consts)
tide_name = np.array([list('M2  ')])
tide_period = np.array([12.420599937439])

# -------------------------------------------------------------------------
g = 9.81
hamp0 = 2.0
cff = np.sqrt(g*h)
cmax0 = hamp0*cff/h
# set variables
hamp = np.ones((consts_num, eta, xi))*hamp0
hpha = np.ones((consts_num, eta, xi))*0.0
cmax = np.ones((consts_num, eta, xi))*cmax0
cmin = np.ones((consts_num, eta, xi))*0.0
cang = np.ones((consts_num, eta, xi))*0.0
cpha = np.ones((consts_num, eta, xi))*0.0

# -------------------------------------------------------------------------
# write tidal information to nc file
# -------------------------------------------------------------------------
# create nc file
fh = nc.Dataset('/Users/cw686/roms_stuff/tide/fjord_tide.nc', 'w')
fh.createDimension('namelen', 4)
fh.createDimension('tide_period', consts_num)
fh.createDimension('eta_rho', eta)
fh.createDimension('xi_rho', xi)

fh.history = 'Tides from TPXO8'
import time
fh.creation_date = time.strftime('%c')
fh.Type = 'ROMS Tidal Forcing File'
fh.Title = 'Forcing for' + dst_grd.name + 'domain'
fh.grid_file = dst_grd.name
fh.Source = 'Analytical'

name_nc = fh.createVariable('tide_name', 'c', ('tide_period', 'namelen'))

period_nc = fh.createVariable('tide_period', 'd', ('tide_period'), fill_value=-9999)
period_nc.field = 'tide_period, scalar'
period_nc.long_name = 'tidal angular period'
period_nc.units = 'hours'

x_nc = fh.createVariable('x_rho', 'd', ('eta_rho', 'xi_rho'))
x_nc.field = 'x_rho, scalar'
x_nc.long_name = 'x coordinates of RHO-points'
x_nc.units = 'm'

y_nc = fh.createVariable('y_rho', 'd', ('eta_rho', 'xi_rho'))
y_nc.field = 'y_rho, scalar'
y_nc.long_name = 'y coordinates of RHO-points'
y_nc.units = 'm'

msk_nc = fh.createVariable('mask_rho', 'd', ('eta_rho', 'xi_rho'))
msk_nc.long_name = 'mask on RHO-points'
msk_nc.option_0 = 'land'
msk_nc.option_1 = 'water'

Eamp_nc = fh.createVariable('tide_Eamp', 'd', ('tide_period', 'eta_rho', 'xi_rho'), fill_value=-9999)
Eamp_nc.field = 'tide_Eamp, scalar'
Eamp_nc.long_name = 'tidal elevation amplitude'
Eamp_nc.units = 'meter'

Ephase_nc = fh.createVariable('tide_Ephase', 'd', ('tide_period', 'eta_rho', 'xi_rho'), fill_value=-9999)
Ephase_nc.field = 'tide_Ephase, scalar'
Ephase_nc.long_name = 'tidal elevation phase angle'
Ephase_nc.units = 'degrees, time of maximum elevation with respect chosen time orgin'

Cmax_nc = fh.createVariable('tide_Cmax', 'd', ('tide_period', 'eta_rho', 'xi_rho'), fill_value=-9999)
Cmax_nc.field = 'tide_Cmax, scalar'
Cmax_nc.long_name = 'maximum tidal current, ellipse semi-major axis'
Cmax_nc.units = 'meter second-1'

Cmin_nc = fh.createVariable('tide_Cmin', 'd', ('tide_period', 'eta_rho', 'xi_rho'), fill_value=-9999)
Cmin_nc.field = 'tide_Cmin, scalar'
Cmin_nc.long_name = 'minimum tidal current, ellipse semi-minor axis'
Cmin_nc.units = 'meter second-1'

Cangle_nc = fh.createVariable('tide_Cangle', 'd', ('tide_period', 'eta_rho', 'xi_rho'), fill_value=-9999)
Cangle_nc.field = 'tide_Cangle, scalar'
Cangle_nc.long_name = 'tidal current inclination angle'
Cangle_nc.units = 'degrees between semi-major axis and East'

Cphase_nc = fh.createVariable('tide_Cphase', 'd', ('tide_period', 'eta_rho', 'xi_rho'), fill_value=-9999)
Cphase_nc.field = 'tide_Cphase, scalar'
Cphase_nc.long_name = 'tidal current phase angle'
Cphase_nc.units = 'degrees, time of maximum velocity'

# -------------------------------------------------------------------------
# write data
name_nc[:, :] = tide_name
period_nc[:] = tide_period

x_nc[:, :] = x
y_nc[:, :] = y
msk_nc[:, :] = msk
Eamp_nc[:, :, :] = hamp
Ephase_nc[:, :, :] = hpha
Cmax_nc[:, :, :] = cmax
Cmin_nc[:, :, :] = cmin
Cangle_nc[:, :, :] = cang
Cphase_nc[:, :, :] = cpha

fh.close()
