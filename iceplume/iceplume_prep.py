import numpy as np
import matplotlib.pyplot as plt
# import netCDF4 as nc
# from datetime import datetime

import pyroms

# get CTD data
from ocean_toolbox import ctd
info = {'data_dir': '/Users/CnWang/Documents/gb_roms/ctd_raw/',
        'file_dir': '/Users/CnWang/Documents/gb_roms/',
        'file_name': 'ctd.nc',
        'sl': 'l',
        'var': ['salt', 'temp', 'o2', 'rho', 'pre', 'fluor', 'tur', 'par'],
        'clim_station': [12],
        'clim_deep_interp': 'yes',
        'filter': 'no',
       }
c = ctd.ctd(info)
c()

stn = c.climatology['station'] == 12
month = 7
depth = -c.climatology['z']
salt = c.climatology['salt'][:, month-1, stn].squeeze()
temp = c.climatology['temp'][:, month-1, stn].squeeze()

# get ROMS grid info
grd = pyroms.grid.get_ROMS_grid('fjord')
zw = grd.vgrid.z_w[:][:, 1, 101]
saltr = np.interp(zw, depth[::-1], salt[::-1])
tempr = np.interp(zw, depth[::-1], temp[::-1])

zw = np.linspace(-40., 0., 41)

# write to output file
fh = open('./iceplume_test_input.txt', 'w')
for i in range(len(zw)):
    fh.write('%f, %f, %f\n' % (zw[i], saltr[i], temp[i]))
fh.close()
