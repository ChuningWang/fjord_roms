from datetime import datetime
import numpy as np
import netCDF4 as nc
import pyroms

# ------------ functionals ---------------------------------------------

# ------------ basic grid parameters -----------------------------------
# grid dimension
Xfjord = 100
Xocean = 50
Yfjord = 500
Yocean = 100

Xpsi = Xfjord + 2*Xocean  # horizontal total
Ypsi = Yfjord + Yocean  # vertical total

Xrho = Xpsi+2
Yrho = Ypsi+2
Xvert = Xpsi+3
Yvert = Ypsi+3

# # total grid length
# Xsize = 100.e3  # m
# Ysize = 150.e3  # m
# 
# # average grid size
# dxa = Xsize/Xpsi  # m
# dya = Ysize/Ypsi  # m

# location of fjord mouth
Ym = 100.e3  # m

# location of sill
Ys = 90.e3  # m

# sill width
Ysw = 5.e3  # m

# location of channel deepening
Yd1 = 60.e3  # m
Yd2 = 70.e3  # m

# fjord width
Xw = 20.e3  # m

# grid spacing parameter
xexp = 0.1
yexp = 0.04

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
grd_name = 'fjord_ana'

# ------------ horizontal grid construction ----------------------------
# X direction
xfjord = np.linspace(-Xw/2., Xw/2., Xfjord+1)
dxfjord = np.diff(xfjord).mean()

dxocean = np.exp(xexp*np.arange(Xocean))*dxfjord
xr = xfjord[-1] + np.cumsum(dxocean)
xl = -xr[::-1]

xvert = np.concatenate((xl, xfjord, xr))
xvert = np.concatenate((np.array([2*xvert[0]-xvert[1]]),
                        xvert,
                        np.array([2*xvert[-1]-xvert[-2]])))

# Y direction
yfjord = np.linspace(0, Ym, Yfjord+1)
dyfjord = np.diff(yfjord).mean()

dyocean = np.exp(yexp*np.arange(Yocean))*dyfjord
yocean = yfjord[-1] + np.cumsum(dyocean)

yvert = np.concatenate((yfjord, yocean))
yvert = np.concatenate((np.array([2*yvert[0]-yvert[1]]),
                        yvert,
                        np.array([2*yvert[-1]-yvert[-2]])))

xxvert, yyvert = np.meshgrid(xvert, yvert)

# generate mask
msk = ((xxvert < -Xw/2.) | (xxvert > Xw/2.)) & (yyvert < Ym)
xxvert = np.ma.masked_where(msk, xxvert)
yyvert = np.ma.masked_where(msk, yyvert)

# ------------ write hgrid ---------------------------------------------
hgrd = pyroms.hgrid.CGrid(xxvert, yyvert)
# Coriolis Parameter
hgrd.f = f0 + hgrd.y_rho*beta

# ------------ vertical grid construction ------------------------------
h = np.zeros(yvert.shape)
# deep water
h[yvert <= Yd1] = Dm 
# shallow water
h[(yvert >= Yd2) & (yvert <= Ym)] = Di
# in-fjord slope
msks = (yvert >= Yd1) & (yvert <= Yd2)
ys = yvert[msks]
h[msks] = Dm - (Dm - Di)*0.5*(1 - np.cos((ys - ys[0])/(ys[-1] - ys[0])*np.pi))
# out-fjord slope
msko = yvert >= Ym
h[msko] = Di + sig*(yvert[msko]-yvert[msko][0])
# sill
msks = (yvert >= Ys-Ysw) * (yvert <= Ys+Ysw)
h[msks] = Di + (Ds - Di)*np.exp(-(yvert[msks]-Ys)**2/(0.1*Ysw**2))

# repmat
h = np.tile(0.5*(h[1:] + h[:-1]), (Xrho, 1)).T
h[hgrd.mask_rho == 0] = hmin

# constrain maximum depth
h[h > 4000] = 4000

# ------------ write vgrid ---------------------------------------------
vgrd = pyroms.vgrid.s_coordinate_4(h, theta_b, theta_s, Tcline, N, hraw=h)

# ------------ write grid ----------------------------------------------
grd = pyroms.grid.ROMS_Grid(grd_name, hgrd, vgrd)
pyroms.grid.write_ROMS_grid(grd, filename='fjord_grd.nc')

