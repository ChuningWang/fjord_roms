from datetime import datetime
import numpy as np
import netCDF4 as nc
import pyroms

# ------------ functionals ---------------------------------------------

# ------------ basic grid parameters -----------------------------------
# grid dimension
Xfjord = 500
Xocean = 50
Yfjord = 100
Yocean = 50

Xpsi = Xfjord + Xocean  # vertical total
Ypsi = Yfjord + 2*Yocean  # horizontal total

Xrho = Xpsi+2
Yrho = Ypsi+2
Xvert = Xpsi+3
Yvert = Ypsi+3

# location of fjord mouth
Xm = 100.e3  # m

# location of sill
Xs = 90.e3  # m

# sill width
Xsw = 5.e3  # m

# location of channel deepening
Xd1 = 60.e3  # m
Xd2 = 70.e3  # m

# fjord width
Yw = 20.e3  # m

# grid spacing parameter
# xinc = 20.0
# yinc = 20.0
xinc = 20.0
yinc = 20.0
xexp = 0.04
yexp = 0.1

# fjord max depth
Dm = 400.

# sill depth
Ds = 30.

# fjord shallow depth
Di = 100.

# Coriolis Parameter
f0 = 1.e-4
beta = 0.

# shelf slope
sig = 0.01

# vertical grid specs
theta_b = 2.0
theta_s = 8.0
Tcline = 10
N = 40
hmin = 10

# grid name
grd_name = 'fjord'

# ------------ horizontal grid construction ----------------------------
# X direction
xfjord = np.linspace(0, Xm, Xfjord+1)
dxfjord = np.diff(xfjord).mean()

# dxocean = np.exp(xexp*np.arange(Xocean))*dxfjord
dxocean = xinc*np.arange(Xocean) + dxfjord
xocean = xfjord[-1] + np.cumsum(dxocean)

xvert = np.concatenate((xfjord, xocean))
xvert = np.concatenate((np.array([2*xvert[0]-xvert[1]]),
                        xvert,
                        np.array([2*xvert[-1]-xvert[-2]])))

# Y direction
yfjord = np.linspace(-Yw/2., Yw/2., Yfjord+1)
dyfjord = np.diff(yfjord).mean()

# dyocean = np.exp(yexp*np.arange(Yocean))*dyfjord
dyocean = yinc*np.arange(Yocean) + dyfjord
yr = yfjord[-1] + np.cumsum(dyocean)
yl = -yr[::-1]

yvert = np.concatenate((yl, yfjord, yr))
yvert = np.concatenate((np.array([2*yvert[0]-yvert[1]]),
                        yvert,
                        np.array([2*yvert[-1]-yvert[-2]])))

# meshgrid it
xxvert, yyvert = np.meshgrid(xvert, yvert)

# generate land mask
msk = np.ones((Yrho, Xrho))
msk0 = ((yyvert < -Yw/2.) | (yyvert > Yw/2.)) & (xxvert < Xm)
msk0 = msk0[:-1,:-1] | msk0[1:,:-1] | msk0[:-1,1:] | msk0[1:,1:]
msk[msk0] = 0
# west boundary
msk[:, 0] = 0

# ------------ write hgrid ---------------------------------------------
hgrd = pyroms.hgrid.CGrid(xxvert, yyvert)
# Coriolis Parameter
hgrd.f = f0 + hgrd.y_rho*beta
hgrd.mask_rho = msk

# ------------ vertical grid construction ------------------------------
h = np.zeros(xvert.shape)
# deep water
h[xvert <= Xd1] = Dm
# shallow water
h[(xvert >= Xd2) & (xvert <= Xm)] = Di
# in-fjord slope
msks = (xvert >= Xd1) & (xvert <= Xd2)
xs = xvert[msks]
h[msks] = Dm - (Dm - Di)*0.5*(1 - np.cos((xs - xs[0])/(xs[-1] - xs[0])*np.pi))
# out-fjord slope
msko = xvert >= Xm
h[msko] = Di + sig*(xvert[msko]-xvert[msko][0])
# sill
msks = (xvert >= Xs-Xsw) * (xvert <= Xs+Xsw)
h[msks] = Di + (Ds - Di)*np.exp(-(xvert[msks]-Xs)**2/(0.1*Xsw**2))

# repmat
h = np.tile(0.5*(h[1:] + h[:-1]), (Yrho, 1))
h[hgrd.mask_rho == 0] = hmin

# constrain maximum depth
h[h > 400] = 400

# ------------ write vgrid ---------------------------------------------
vgrd = pyroms.vgrid.s_coordinate_4(h, theta_b, theta_s, Tcline, N, hraw=h)

# ------------ write grid ----------------------------------------------
grd = pyroms.grid.ROMS_Grid(grd_name, hgrd, vgrd)
pyroms.grid.write_ROMS_grid(grd, filename='fjord_grd.nc')

