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

eta, xi = x.shape
fnum = 2
snum = 24*2
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
    fin = nc.Dataset('./data/run_iceplume/fjord_his_%04d.nc' % (i+1))
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

    fin = nc.Dataset('./data/run_noiceplume/fjord_his_%04d.nc' % (i+1))
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

# x = x[1:-1, 1:-1]
# y = y[1:-1, 1:-1]
# zr = zr[:, 1:-1, 1:-1]
# zw = zw[:, 1:-1, 1:-1]
# 
# zeta = zeta[:, 1:-1, 1:-1]
# salt = salt[:, :, 1:-1, 1:-1]
# temp = temp[:, :, 1:-1, 1:-1]
# u = 0.5*(u[:, :, 1:-1, 1:] + u[:, :, 1:-1, :-1])
# v = 0.5*(v[:, :, 1:, 1:-1] + v[:, :, :-1, 1:-1])
# w = 0.5*(w[:, 1:, 1:-1, 1:-1] + w[:, :-1, 1:-1, 1:-1])
# 
# zeta2 = zeta2[:, 1:-1, 1:-1]
# salt2 = salt2[:, :, 1:-1, 1:-1]
# temp2 = temp2[:, :, 1:-1, 1:-1]
# u2 = 0.5*(u2[:, :, 1:-1, 1:] + u2[:, :, 1:-1, :-1])
# v2 = 0.5*(v2[:, :, 1:, 1:-1] + v2[:, :, :-1, 1:-1])
# w2 = 0.5*(w2[:, 1:, 1:-1, 1:-1] + w2[:, :-1, 1:-1, 1:-1])

w = 0.25*(w[:, 1:, :, 1:] + w[:, 1:, :, :-1] +
          w[:, :-1, :, 1:] + w[:, :-1, :, :-1])
w2 = 0.25*(w2[:, 1:, :, 1:] + w2[:, 1:, :, :-1] +
           w2[:, :-1, :, 1:] + w2[:, :-1, :, :-1])

# Hovmuller diagram
fig, axs = plt.subplots(2)
for i in range(nt):
    # pc1 = axs[0].pcolor(x[3, :], zr[:, 3, 0], salt[i, :, 3, :],
    #                     vmin=15, vmax=30, cmap = cm.haline)
    # pc2 = axs[1].pcolor(x[3, :], zr[:, 3, 0], salt2[i, :, 3, :],
    #                     vmin=15, vmax=30, cmap = cm.haline)
    pc1 = axs[0].pcolor(x[3, :], zr[:, 3, 0], dye[i, :, 3, :],
                        vmin=0, vmax=1, cmap = cm.matter)
    pc2 = axs[1].pcolor(x[3, :], zr[:, 3, 0], dye2[i, :, 3, :],
                        vmin=0, vmax=1, cmap = cm.matter)

    qv1 = axs[0].quiver(xu[3, :], zr[:, 3, 0],
                        u[i, :, 3, :], w[i, :, 3, :],
                        scale=uscale)
    qv2 = axs[1].quiver(xu[3, :], zr[:, 3, 0],
                        u2[i, :, 3, :], w2[i, :, 3, :],
                        scale=uscale)
    qvkey = plt.quiverkey(qv1, 0.9, 0.1, 0.5, r'0.5 ms$^{-1}$')

    cbar_ax = fig.add_axes([0.92, 0.11, 0.02, 0.77])
    cb = fig.colorbar(pc1, cax=cbar_ax)

    plt.savefig('./figs/dye_%04d.png' % (i))
    pc1.remove()
    pc2.remove()
    qv1.remove()
    qv2.remove()
    qvkey.remove()
    fig.delaxes(cbar_ax)

    pc1 = axs[0].pcolor(x[3, :], zr[:, 3, 0], salt[i, :, 3, :],
                        vmin=15, vmax=30, cmap = cm.haline)
    pc2 = axs[1].pcolor(x[3, :], zr[:, 3, 0], salt2[i, :, 3, :],
                        vmin=15, vmax=30, cmap = cm.haline)
    # pc1 = axs[0].pcolor(x[3, :], zr[:, 3, 0], dye[i, :, 3, :],
    #                     vmin=0, vmax=1, cmap = cm.matter)
    # pc2 = axs[1].pcolor(x[3, :], zr[:, 3, 0], dye2[i, :, 3, :],
    #                     vmin=0, vmax=1, cmap = cm.matter)

    qv1 = axs[0].quiver(xu[3, :], zr[:, 3, 0],
                        u[i, :, 3, :], w[i, :, 3, :],
                        scale=uscale)
    qv2 = axs[1].quiver(xu[3, :], zr[:, 3, 0],
                        u2[i, :, 3, :], w2[i, :, 3, :],
                        scale=uscale)
    qvkey = plt.quiverkey(qv1, 0.9, 0.1, 0.5, r'0.5 ms$^{-1}$')

    cbar_ax = fig.add_axes([0.92, 0.11, 0.02, 0.77])
    cb = fig.colorbar(pc1, cax=cbar_ax)

    plt.savefig('./figs/salt_%04d.png' % (i))
    pc1.remove()
    pc2.remove()
    qv1.remove()
    qv2.remove()
    qvkey.remove()
    fig.delaxes(cbar_ax)

plt.close()
