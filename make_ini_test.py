
from datetime import datetime

import numpy as np
import scipy as sp
from scipy import io
import netCDF4 as nc

import pyroms
import pyroms_toolbox

# ----------------- functionals ---------------------------------------------
def get_z(h, hc, N, s_rho, Cs_r, zeta, Vtrans, zice):

    hwater = h - np.abs(zice)
    z_r = np.empty((N,) + h.shape, 'd')
    if Vtrans == 1:
        for k in range(N):
            z0 = hc * s_rho[k] + (hwater - hc) * Cs_r[k]
            z_r[k, :] = z0 + zeta * (1.0 + z0 / hwater)
            z_r[k, :] = z_r[k, :] - np.abs(zice)
    elif Vtrans == 2 or Vtrans == 4 or Vtrans == 5:
        for  k in range(N):
            z0 = (hc * s_rho[k] + hwater * Cs_r[k]) / (hc + hwater)
            z_r[k, :] = zeta + (zeta + hwater) * z0
            z_r[k, :] = z_r[k, :] - np.abs(zice)
    return z_r

# ------------ read grid -----------------------------------------------
grd = pyroms.grid.get_ROMS_grid('fjord_test')
msk = grd.hgrid.mask
zr = grd.vgrid.z_r[:]
N, eta, xi = zr.shape

fh = nc.Dataset('./fjord_grd_test.nc')
zice = fh.variables['zice'][:]
fh.close()

zr = get_z(grd.vgrid.h,
           grd.vgrid.hc, grd.vgrid.N, grd.vgrid.s_rho, grd.vgrid.Cs_r,
           np.zeros(zice.shape), grd.vgrid.Vtrans,
           zice)

# ------------ generate ic file ----------------------------------------
ic_file = './fjord_ic_test.nc'
class ocean_time_info(object):
    pass
ocean_time = ocean_time_info()
ocean_time.long_name = 'seconds since 00-00-00'
ocean_time.units = 'second'

pyroms_toolbox.nc_create_roms_file(ic_file, grd, ocean_time, geogrid=False)

# ------------ write ic values -----------------------------------------
z_pyc = 10.
zeta = 0.
ubar = 0.
vbar = 0.
u = 0.
v = 0.
temp = np.zeros((1, grd.vgrid.N) + grd.vgrid.h.shape)
salt = np.zeros((1, grd.vgrid.N) + grd.vgrid.h.shape)
temp[0, :] = 4. + 2.*(np.tanh(0.025*np.pi*(zr+z_pyc))+1)
salt[0, :] = 30. - 10.*(np.tanh(0.025*np.pi*(zr+z_pyc))+1)
dye_01 = 0.
spval = -1.0e20

# zeta = 0
# ubar = 0
# vbar = 0
# u = 0
# v = 0
# dye_01 = 0
# 
# data = io.loadmat('../data/carroll_2017_JGR_oceans_initialTS.mat')
# depth_raw = np.array(data['depth']).squeeze()
# temp_raw = np.array(data['temperature']).squeeze()
# salt_raw = np.array(data['salinity']).squeeze()
# depth_raw = depth_raw[::-1]
# temp_raw = temp_raw[::-1]
# salt_raw = salt_raw[::-1]
# 
# temp = np.zeros((1, N, eta, xi))
# salt = np.zeros((1, N, eta, xi))
# for ieta in range(eta):
#     for ixi in range(xi):
#         temp[0, :, ieta, ixi] = np.interp(zr[:, ieta, ixi], depth_raw, temp_raw)
#         salt[0, :, ieta, ixi] = np.interp(zr[:, ieta, ixi], depth_raw, salt_raw)
# 
# spval = -1.0e20
# ------------ write to ic file ----------------------------------------
fh = nc.Dataset(ic_file, 'r+')

fh.variables['ocean_time'][:] = 0

fh.createVariable('zeta', 'f8', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=spval)
fh.variables['zeta'].long_name = 'free-surface'
fh.variables['zeta'].units = 'meter'
fh.variables['zeta'].field = 'free-surface, scalar, series'
fh.variables['zeta'][:] = zeta

fh.createVariable('temp', 'f8', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=spval)
fh.variables['temp'].long_name = 'potential temperature'
fh.variables['temp'].units = 'Celsius'
fh.variables['temp'].field = 'temperature, scalar, series'
fh.variables['temp'][:] = temp

fh.createVariable('salt', 'f8', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=spval)
fh.variables['salt'].long_name = 'salinity'
fh.variables['salt'].units = 'PSU'
fh.variables['salt'].field = 'salinity, scalar, series'
fh.variables['salt'][:] = salt

fh.createVariable('dye_01', 'f8', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=spval)
fh.variables['dye_01'].long_name = 'river dye concentration'
fh.variables['dye_01'].units = ' '
fh.variables['dye_01'].field = 'dye_01, scalar, series'
fh.variables['dye_01'][:] = dye_01

fh.createVariable('u', 'f8', ('ocean_time', 's_rho', 'eta_u', 'xi_u'), fill_value=spval)
fh.variables['u'].long_name = '3D u-momentum component'
fh.variables['u'].units = 'meter second-1'
fh.variables['u'].field = 'u-velocity, scalar, series'
fh.variables['u'][:] = u

fh.createVariable('v', 'f8', ('ocean_time', 's_rho', 'eta_v', 'xi_v'), fill_value=spval)
fh.variables['v'].long_name = '3D v-momentum component'
fh.variables['v'].units = 'meter second-1'
fh.variables['v'].field = 'v-velocity, scalar, series'
fh.variables['v'][:] = v

fh.createVariable('ubar', 'f8', ('ocean_time', 'eta_u', 'xi_u'), fill_value=spval)
fh.variables['ubar'].long_name = '2D u-momentum component'
fh.variables['ubar'].units = 'meter second-1'
fh.variables['ubar'].field = 'ubar-velocity, scalar, series'
fh.variables['ubar'][:] = ubar

fh.createVariable('vbar', 'f8', ('ocean_time', 'eta_v', 'xi_v'), fill_value=spval)
fh.variables['vbar'].long_name = '2D v-momentum component'
fh.variables['vbar'].units = 'meter second-1'
fh.variables['vbar'].field = 'vbar-velocity, scalar, series'
fh.variables['vbar'][:] = vbar

fh.close()
