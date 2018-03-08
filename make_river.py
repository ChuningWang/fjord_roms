import numpy as np
import netCDF4 as nc
from datetime import datetime, timedelta

import pyroms

# ---------------- load grid ----------------------------------
grd1 = 'fjord'
grd = pyroms.grid.get_ROMS_grid(grd1)
msk = grd.hgrid.mask[:, 1]

xpos = np.where(msk == 1)[0]
ypos = np.ones(xpos.shape)
nr = len(ypos)
rdir = np.zeros(nr)
river = np.ones(nr)
rtime = np.array([0., 366.])
rtrs = np.ones((2, nr))*50.
rtemp = np.ones((2, grd.vgrid.N, nr))*10.
rsalt = np.ones((2, grd.vgrid.N, nr))*0.
rdye1 = np.ones((2, grd.vgrid.N, nr))*1.

# v_shape = np.zeros((grd.vgrid.N, nr))
# v_shape[-1, :] = 1

v_shape = np.ones((grd.vgrid.N, nr))
v_shape = v_shape/grd.vgrid.N

# create file with all the objects
fout = nc.Dataset('fjord_river.nc', 'w')
fout.type = 'ROMS RIVERS file'
fout.title = 'Fjord test'
fout.source = 'Analytical'

fout.createDimension('river_time', None)
fout.createDimension('river', nr)
fout.createDimension('s_rho', grd.vgrid.N)

times = fout.createVariable('river_time', 'f8', ('river_time'))
times.units = 'days since 1900-01-01 00:00:00'
times.long_name = 'river runoff time'
fout.variables['river_time'][:] = rtime

rivers = fout.createVariable('river', 'i4', ('river'))
rivers.long_name = 'river runoff identification number'
fout.variables['river'][:] = river

eta = fout.createVariable('river_Eposition', 'i4', ('river'))
eta.long_name = 'river ETA-position at RHO-points'
fout.variables['river_Eposition'][:] = xpos

xi = fout.createVariable('river_Xposition', 'i4', ('river'))
xi.long_name = 'river XI-position at RHO-points'
fout.variables['river_Xposition'][:] = ypos

dirs = fout.createVariable('river_direction', 'i4', ('river'))
dirs.long_name = 'river runoff direction'
fout.variables['river_direction'][:] = rdir

trans = fout.createVariable('river_transport', 'f8', ('river_time', 'river'))
trans.long_name = 'river runoff vertically integrated mass transport'
trans.units = 'meter3 second-1'
trans.time = 'river_time'
fout.variables['river_transport'][:] = rtrs

vshape = fout.createVariable('river_Vshape', 'f8', ('s_rho', 'river'))
vshape.long_name = 'river runoff mass transport vertical profile'
fout.variables['river_Vshape'][:] = v_shape

temp = fout.createVariable('river_temp', 'f8', ('river_time', 's_rho', 'river'))
temp.long_name = 'river runoff potential temperature'
temp.units = 'Celsius'
temp.time = 'river_time'
fout.variables['river_temp'][:] = rtemp

salt = fout.createVariable('river_salt', 'f8', ('river_time', 's_rho', 'river'))
salt.long_name = 'river runoff salinity'
salt.time = 'river_time'
fout.variables['river_salt'][:] = rsalt

dye1 = fout.createVariable('river_dye_01', 'f8', ('river_time', 's_rho', 'river'))
dye1.long_name = 'river runoff dye 01 concentration'
dye1.time = 'river_time'
fout.variables['river_dye_01'][:] = rdye1
fout.close()
