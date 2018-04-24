import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

from cmocean import cm
import pyroms

uscale = 5
# ----------------- load grid -----------------------------------------------
grd = pyroms.grid.get_ROMS_grid('fjord_test')
x = grd.hgrid.x_rho
y = grd.hgrid.y_rho
xvert = grd.hgrid.x_vert
yvert = grd.hgrid.y_vert
xu = grd.hgrid.x_u
yu = grd.hgrid.y_u
xv = grd.hgrid.x_v
yv = grd.hgrid.y_v
zr = grd.vgrid.z_r[:]
zw = grd.vgrid.z_w[:]
msk = grd.hgrid.mask_rho
msku = grd.hgrid.mask_u
mskv = grd.hgrid.mask_v
N = grd.vgrid.N

jslice = 3

xvert_t = np.tile(xvert[jslice, :], (N+1, 1))
zw_t = np.zeros(xvert_t.shape)
zw_t[:, 1:-1] = 0.5 * (zw[:, jslice, 1:] + zw[:, jslice, :-1])
zw_t[:, 1] = zw[:, jslice, 1]
zw_t[:, -1] = zw[:, jslice, -1]

xu_t = np.tile(xu[jslice, :], (N, 1))
zr_t = 0.5 * (zr[:, jslice, 1:] + zr[:, jslice, :-1])

eta, xi = x.shape
fnum = 10
snum = 24*1
nt = fnum*snum+1
time = np.zeros((nt))

zeta = np.zeros((nt, eta, xi))
salt = np.zeros((nt, N, eta, xi))
temp = np.zeros((nt, N, eta, xi))
dye = np.zeros((nt, N, eta, xi))
u = np.zeros((nt, N, eta, xi-1))
v = np.zeros((nt, N, eta-1, xi))
w = np.zeros((nt, N+1, eta, xi))

zeta2 = np.zeros((nt, eta, xi))
salt2 = np.zeros((nt, N, eta, xi))
temp2 = np.zeros((nt, N, eta, xi))
dye2 = np.zeros((nt, N, eta, xi))
u2 = np.zeros((nt, N, eta, xi-1))
v2 = np.zeros((nt, N, eta-1, xi))
w2 = np.zeros((nt, N+1, eta, xi))

msk3d = np.tile(msk == 0, (nt, 1, 1))
msk4d = np.tile(msk == 0, (nt, N, 1, 1))
mskw4d = np.tile(msk == 0, (nt, N+1, 1, 1))
msku3d = np.tile(msku == 0, (nt, 1, 1))
msku4d = np.tile(msku == 0, (nt, N, 1, 1))
mskv3d = np.tile(mskv == 0, (nt, 1, 1))
mskv4d = np.tile(mskv == 0, (nt, N, 1, 1))

for i in range(fnum):
    fin = nc.Dataset('/Users/cw686/roms_stuff/outputs/iceplume/fjord_his_%04d.nc' % (i+1))
    if i == 0:
        time[:snum+1] = fin.variables['ocean_time'][:]
        zeta[:snum+1, :, :] = fin.variables['zeta'][:]
        salt[:snum+1, :, :, :] = fin.variables['salt'][:]
        temp[:snum+1, :, :, :] = fin.variables['temp'][:]
        dye[:snum+1, :, :, :] = fin.variables['dye_01'][:]
        u[:snum+1, :, :, :] = fin.variables['u'][:]
        v[:snum+1, :, :, :] = fin.variables['v'][:]
        w[:snum+1, :, :, :] = fin.variables['w'][:]
    else:
        time[i*snum+1:(i+1)*snum+1] = fin.variables['ocean_time'][:]
        zeta[i*snum+1:(i+1)*snum+1, :, :] = fin.variables['zeta'][:]
        salt[i*snum+1:(i+1)*snum+1, :, :, :] = fin.variables['salt'][:]
        temp[i*snum+1:(i+1)*snum+1, :, :, :] = fin.variables['temp'][:]
        dye[i*snum+1:(i+1)*snum+1, :, :, :] = fin.variables['dye_01'][:]
        u[i*snum+1:(i+1)*snum+1, :, :, :] = fin.variables['u'][:]
        v[i*snum+1:(i+1)*snum+1, :, :, :] = fin.variables['v'][:]
        w[i*snum+1:(i+1)*snum+1, :, :, :] = fin.variables['w'][:]
    fin.close()

    fin = nc.Dataset('/Users/cw686/roms_stuff/outputs/ref/fjord_his_%04d.nc' % (i+1))
    if i == 0:
        zeta2[:snum+1, :, :] = fin.variables['zeta'][:]
        salt2[:snum+1, :, :, :] = fin.variables['salt'][:]
        temp2[:snum+1, :, :, :] = fin.variables['temp'][:]
        dye2[:snum+1, :, :, :] = fin.variables['dye_01'][:]
        u2[:snum+1, :, :, :] = fin.variables['u'][:]
        v2[:snum+1, :, :, :] = fin.variables['v'][:]
        w2[:snum+1, :, :, :] = fin.variables['w'][:]
    else:
        zeta2[i*snum+1:(i+1)*snum+1, :, :] = fin.variables['zeta'][:]
        salt2[i*snum+1:(i+1)*snum+1, :, :, :] = fin.variables['salt'][:]
        temp2[i*snum+1:(i+1)*snum+1, :, :, :] = fin.variables['temp'][:]
        dye2[i*snum+1:(i+1)*snum+1, :, :, :] = fin.variables['dye_01'][:]
        u2[i*snum+1:(i+1)*snum+1, :, :, :] = fin.variables['u'][:]
        v2[i*snum+1:(i+1)*snum+1, :, :, :] = fin.variables['v'][:]
        w2[i*snum+1:(i+1)*snum+1, :, :, :] = fin.variables['w'][:]
    fin.close()

time = time/3600.

zeta = np.ma.masked_where(msk3d, zeta)
temp = np.ma.masked_where(msk4d, temp)
salt = np.ma.masked_where(msk4d, salt)
dye = np.ma.masked_where(msk4d, dye)
u = np.ma.masked_where(msku4d, u)
v = np.ma.masked_where(mskv4d, v)
# w = np.ma.masked_where(mskw4d, w)

zeta2 = np.ma.masked_where(msk3d, zeta2)
temp2 = np.ma.masked_where(msk4d, temp2)
salt2 = np.ma.masked_where(msk4d, salt2)
dye2 = np.ma.masked_where(msk4d, dye2)
u2 = np.ma.masked_where(msku4d, u2)
v2 = np.ma.masked_where(mskv4d, v2)
# w2 = np.ma.masked_where(mskw4d, w2)

w = 0.25*(w[:, 1:, :, 1:] + w[:, 1:, :, :-1] +
          w[:, :-1, :, 1:] + w[:, :-1, :, :-1])
w2 = 0.25*(w2[:, 1:, :, 1:] + w2[:, 1:, :, :-1] +
           w2[:, :-1, :, 1:] + w2[:, :-1, :, :-1])

# Along track transects
fig, axs = plt.subplots(4)
fig.subplots_adjust(hspace=0.05, right=0.85)
axs[3].set_xlabel('X [m]')
axs[3].set_ylabel('Z [m]')
# axs[0].set_ylabel('Z [m]')
axs[0].text(8000, -25, 'ICEPLUME')
axs[1].text(8000, -25, 'NO ICEPLUME')
axs[2].text(8000, -25, 'ICEPLUME')
axs[3].text(8000, -25, 'NO ICEPLUME')
axs[0].set_xticklabels([''])
axs[1].set_xticklabels([''])
axs[2].set_xticklabels([''])
axs[0].set_yticks([-50, -100, -150, -200])
axs[0].set_yticklabels(['50', '100', '150', '200'])
axs[1].set_yticks([-50, -100, -150, -200])
axs[1].set_yticklabels(['50', '100', '150', '200'])
axs[2].set_yticks([-50, -100, -150, -200])
axs[2].set_yticklabels(['50', '100', '150', '200'])
axs[3].set_yticks([-50, -100, -150, -200])
axs[3].set_yticklabels(['50', '100', '150', '200'])
for i in range(nt):
    pc1 = axs[0].pcolor(xvert_t, zw_t, dye[i, :, 3, :],
                        vmin=0, vmax=.5, cmap = cm.matter)
    pc2 = axs[1].pcolor(xvert_t, zw_t, dye2[i, :, 3, :],
                        vmin=0, vmax=.5, cmap = cm.matter)

    qv1 = axs[0].quiver(xu_t, zr_t,
                        u[i, :, 3, :], w[i, :, 3, :],
                        scale=uscale)
    qv2 = axs[1].quiver(xu_t, zr_t,
                        u2[i, :, 3, :], w2[i, :, 3, :],
                        scale=uscale)
    # qvkey = plt.quiverkey(qv1, 0.9, 0.1, 0.5, r'0.5 ms$^{-1}$')

    pc3 = axs[2].pcolor(xvert_t, zw_t, salt[i, :, 3, :],
                        vmin=10, vmax=30, cmap = cm.haline)
    pc4 = axs[3].pcolor(xvert_t, zw_t, salt2[i, :, 3, :],
                        vmin=10, vmax=30, cmap = cm.haline)

    qv3 = axs[2].quiver(xu_t, zr_t,
                        u[i, :, 3, :], w[i, :, 3, :],
                        scale=uscale)
    qv4 = axs[3].quiver(xu_t, zr_t,
                        u2[i, :, 3, :], w2[i, :, 3, :],
                        scale=uscale)

    qvkey = plt.quiverkey(qv2, 0.9, 0.1, 0.5, r'0.5 ms$^{-1}$')

    if i==0:
        cbar_ax1 = fig.add_axes([0.87, 0.51, 0.02, 0.37])
        cb1 = fig.colorbar(pc1, cax=cbar_ax1,
                           ticks=[0, 0.1, 0.2, 0.3, 0.4, 0.5])
        cb1.set_label('Dye')

        cbar_ax2 = fig.add_axes([0.87, 0.11, 0.02, 0.37])
        cb2 = fig.colorbar(pc3, cax=cbar_ax2,
                           ticks=[10, 15, 20, 25, 30])
        cb2.set_label('Salinity [PSU]')

    ttl = axs[0].set_title(grd.name + '_Hour_' + str(time[i]))

    plt.savefig('./figs/his_dye_salt_%04d.png' % (i))
    pc1.remove()
    pc2.remove()
    qv1.remove()
    qv2.remove()
    qvkey.remove()
    # fig.delaxes(cbar_ax)

    # plt.savefig('./figs/salt_%04d.png' % (i))
    pc3.remove()
    pc4.remove()
    qv3.remove()
    qv4.remove()
    # qvkey.remove()
    # fig.delaxes(cbar_ax)

plt.close()
