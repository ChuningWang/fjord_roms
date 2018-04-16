from datetime import datetime
import numpy as np
import netCDF4 as nc
import pyroms

# ------------ functionals ---------------------------------------------

# ------------ basic grid parameters -----------------------------------
# horizontal resolution
dx = 200.
dy = 200.

# channel dimensions
L = 60000
W = 20000

# total dimensions
Lt = 125000
Wt = 160000

# ocean dimensions
Locean = Lt - L
Wocean = 0.5 * (Wt - W)

# grid dimension
Xfjord = int(L/dx)
Xocean = int(360-Xfjord)
Yfjord = int(W/dy)
Yocean = int((270-Yfjord)/2)

Xpsi = Xfjord + Xocean  # vertical total
Ypsi = Yfjord + 2*Yocean  # horizontal total

Xrho = Xpsi+2
Yrho = Ypsi+2
Xvert = Xpsi+3
Yvert = Ypsi+3

# sill width
Xsw = 5.e3  # m

# location of sill
Xs = L-Xsw  # m

# location of channel deepening
Xd1 = 30.e3  # m
Xd2 = 40.e3  # m

# # grid spacing parameter
# xinc = 20.0
# yinc = 20.0
# xexp = 0.04
# yexp = 0.1

# fjord max depth
Dm = 800.

# sill depth
Ds = 400.

# fjord shallow depth
Di = 800.

# Coriolis Parameter
f0 = 1.e-4
beta = 0.

# shelf slope
sig = 0.00

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
xfjord = np.linspace(0, L, Xfjord+1)
dxfjord = dx

# dxocean = np.exp(xexp*np.arange(Xocean))*dxfjord
# dxocean = xinc*np.arange(Xocean) + dxfjord

ddxocean = (2*Locean/Xocean-2*dx)/(1+Xocean)
dxocean = dx + (np.arange(Xocean) + 1)*ddxocean

xocean = xfjord[-1] + np.cumsum(dxocean)

xvert = np.concatenate((xfjord, xocean))
xvert = np.concatenate((np.array([2*xvert[0]-xvert[1]]),
                        xvert,
                        np.array([2*xvert[-1]-xvert[-2]])))

# Y direction
yfjord = np.linspace(-W/2., W/2., Yfjord+1)
dyfjord = dy

# dyocean = np.exp(yexp*np.arange(Yocean))*dyfjord
# dyocean = yinc*np.arange(Yocean) + dyfjord

ddyocean = (2*Wocean/Yocean-2*dy)/(1+Yocean)
dyocean = dy + (np.arange(Yocean) + 1)*ddyocean

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
msk0 = ((yyvert < -W/2.) | (yyvert > W/2.)) & (xxvert < L)
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
h[(xvert >= Xd2) & (xvert <= L)] = Di
# in-fjord slope
msks = (xvert >= Xd1) & (xvert <= Xd2)
xs = xvert[msks]
h[msks] = Dm - (Dm - Di)*0.5*(1 - np.cos((xs - xs[0])/(xs[-1] - xs[0])*np.pi))
# out-fjord slope
msko = xvert >= L
h[msko] = Di + sig*(xvert[msko]-xvert[msko][0])
# sill
msks = (xvert >= Xs-Xsw) * (xvert <= Xs+Xsw)
h[msks] = Di + (Ds - Di)*np.exp(-(xvert[msks]-Xs)**2/(0.1*Xsw**2))

# repmat
h = np.tile(0.5*(h[1:] + h[:-1]), (Yrho, 1))
h[hgrd.mask_rho == 0] = hmin

# constrain maximum depth
h[h > Dm] = Dm

# ------------ write vgrid ---------------------------------------------
vgrd = pyroms.vgrid.s_coordinate_4(h, theta_b, theta_s, Tcline, N, hraw=h)

# ------------ write grid ----------------------------------------------
grd = pyroms.grid.ROMS_Grid(grd_name, hgrd, vgrd)
pyroms.grid.write_ROMS_grid(grd, filename='fjord_grd.nc')

