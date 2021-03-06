"""
Plot map view of variables in Glacier Bay.
"""

# --------------------- load modules --------------------------------------
import sys
import glob
from datetime import datetime, timedelta

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.basemap import Basemap

import cmocean
import pyroms
import pyroms_toolbox as prt

# --------------------- input arguments -----------------------------------
# my inputs
my_year = 2008
grd1 = 'fjord'
ftype = 'his'
varlist = ['zeta', 'sbar']
depth = 1
dd = 2
plt_uv = 1
plt_contourf = 0
uscale = 10

outputs_dir = '/glade/scratch/chuning/tmpdir_fjord/outputs/'
fig_dir = '/glade/p/work/chuning/fjord/figs/zview/' + ftype + '/'

flist = sorted(glob.glob(outputs_dir + '*' + ftype + '*.nc'))

# dicts for variable clims, colormaps, and other properties
clim = {'temp': [2, 10],
        'salt': [15, 30],
        'zeta': [-3., 3.],
        'tbar': [2, 10],
        'sbar': [25, 30],
        'wetdry_mask_rho': [0, 1],
        'dye_01': [0, 1],
        'dye_02': [0, 1],
        'dye_03': [0, 1],
        'tke': [1e-5, 1e0],
        'gls': [1e-6, 1e-4],
       }

cmap = {'temp': cmocean.cm.thermal,
        'salt': cmocean.cm.haline,
        'zeta': cmocean.cm.balance,
        'tbar': cmocean.cm.thermal,
        'sbar': cmocean.cm.haline,
        'wetdry_mask_rho': cmocean.cm.balance,
        'dye_01': cmocean.cm.matter,
        'dye_02': cmocean.cm.matter,
        'dye_03': cmocean.cm.matter,
        'tke': cmocean.cm.matter,
        'gls': cmocean.cm.matter,
       }

var_2d = ['zeta', 'wetdry_mask_rho', 'sbar', 'tbar']
var_omega = ['tke', 'gls']
var_log = ['tke', 'gls']

lat_min = 57.75
lat_max = 59.25
lat_0 = 0.5 * (lat_min + lat_max)

lon_min = -137.5
lon_max = -135.0
lon_0 = 0.5 * (lon_min + lon_max)

# --------------------- data preparation ----------------------------------
depth = -abs(depth)
grd = pyroms.grid.get_ROMS_grid(grd1)
x = grd.hgrid.x_rho
y = grd.hgrid.y_rho
mask = grd.hgrid.mask_rho
zw = grd.vgrid.z_w[:]
dz = np.diff(zw, axis=0)
h = grd.vgrid.h

# --------------------- make plots ----------------------------------------
# plot options
plt.switch_backend('Agg')
try:
    plt.style.use('classic')
except:
    pass

# draw background
fig, ax = plt.subplots()
ax.set_xlim(0, 120000)
ax.set_ylim(-20000, 20000)

for var in varlist:

    if var in ['sbar']:
        var2 = 'salt'
    elif var in ['tbar']:
        var2 = 'temp'
    else:
        var2 = var

    if var in clim.keys():
        clim_var = clim[var]
        cmap_var = cmap[var]
    else:
        clim_var = [0, 1]
        cmap_var = cmocean.cm.matter

    if var in var_log:
        clevs = np.logspace(clim_var[0], clim_var[1], 3)
    else:
        clevs = np.linspace(clim_var[0], clim_var[1], 5)

    if var in var_2d:
        uvar = 'ubar'
        vvar = 'vbar'
    else:
        uvar = 'u'
        vvar = 'v'

    for fn in flist:

        tag = fn.split('/')[-1].split('.nc')[0]
        print('processing ' + var + ' for ' + fn + ' ...')

        fh = nc.Dataset(fn)
        t = fh.variables['ocean_time'][:]
        tunit = (fh.variables['ocean_time']).units
        data = fh.variables[var2][:]
        if var in ['sbar', 'tbar']:
            data = np.sum(data*dz, axis=1)/np.tile(h, (len(t), 1, 1))

        if plt_uv == 1:
            u = fh.variables[uvar][:]
            v = fh.variables[vvar][:]
        fh.close()

        if var in var_omega:
            data = 0.5*(data[:, :-1, :, :] + data[:, 1:, :, :])

        for tindex in range(len(t)):
            ttag = nc.num2date(t[tindex], tunit).strftime("%Y-%m-%d_%H:%M:%S")

            if var not in var_2d:
                # get z slice
                zslice, _, _ = pyroms.tools.zslice(data[tindex, :, :, :], depth, grd, \
                                          Cpos='rho', vert=True)
                if plt_uv == 1:
                    # get u and v slice at requested depth
                    zsliceu, _, _ = pyroms.tools.zslice(u[tindex, :, :, :], depth, grd, Cpos='u')
                    zslicev, _, _ = pyroms.tools.zslice(v[tindex, :, :, :], depth, grd, Cpos='v')

            else:
                zslice = data[tindex, :, :]
                if plt_uv == 1:
                    zsliceu = u[tindex, :, :]
                    zslicev = v[tindex, :, :]

            zslice = np.ma.masked_where(mask == 0, zslice)
            if plt_uv == 1:
                # average field at rho point position
                zsliceu = 0.5 * (zsliceu[:,:-1] + zsliceu[:,1:])
                zsliceu = zsliceu[:,np.r_[0,:np.size(zsliceu,1),-1]]
                zsliceu = np.ma.masked_where(mask == 0, zsliceu)
                zsliceu = np.ma.masked_where(zsliceu >= 1000, zsliceu)
                zslicev = 0.5 * (zslicev[:-1,:] + zslicev[1:,:])
                zslicev = zslicev[np.r_[0,:np.size(zslicev,0),-1],:]
                zslicev = np.ma.masked_where(mask == 0, zslicev)
                zslicev = np.ma.masked_where(zslicev >= 1000, zslicev)
                U = zsliceu + 1j * zslicev
                # rotate velocity vector according to grid angle
                U = U * np.exp(1j * grd.hgrid.angle_rho)

            if var in var_log:
                pc = ax.pcolormesh(x, y, zslice, norm=LogNorm(vmin=clim_var[0], vmax=clim_var[1]), cmap=cmap_var)
            else:
                pc = ax.pcolormesh(x, y, zslice, cmap=cmap_var)
                pc.set_clim(clim_var[0], clim_var[1])

            if plt_uv == 1:
                qv = ax.quiver(x[::dd,::dd], y[::dd,::dd], \
                              np.real(U[::dd,::dd]), np.imag(U[::dd,::dd]), \
                              scale = uscale, width=0.001)
                qvkey = plt.quiverkey(qv, 0.2, 0.1, 1, r'1 ms$^{-1}$')

            cbar_ax = fig.add_axes([0.9, 0.12, 0.02, 0.76])
            cb = fig.colorbar(pc, cax=cbar_ax)

            if plt_contourf == 1:
                varc = ax.contour(x, y, zslice, clevs, linestyle='--', linewidths=.4, colors='w')
                # varcl = plt.clabel(varc, fontsize=5)

            if var not in var_2d:
                ttl = ax.set_title(var + '_' + str(int(-depth)) + 'm_' + grd.name + '_' + ftype + '_' + ttag)
                fig.savefig(fig_dir + var + '_' + str(int(abs(depth))) + 'm_' + grd.name + '_' + ftype + '_' + ttag + '.png')
            else:
                ttl = ax.set_title(var + '_' + ttag)
                fig.savefig(fig_dir + var + '_' + grd.name + '_' + ftype + '_' + ttag + '.png')

            pc.remove()
            fig.delaxes(cbar_ax)
            # cb.remove()
            if plt_contourf == 1:
                for cc in varc.collections:
                    cc.remove()
            if plt_uv == 1:
                qv.remove()
                qvkey.remove()

plt.close()
