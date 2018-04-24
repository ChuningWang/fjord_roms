from datetime import datetime
import numpy as np
import netCDF4 as nc
import pyroms

# ------------ functionals ---------------------------------------------

# ------------ basic grid parameters -----------------------------------
# grid dimension
Xfjord = 100
Yfjord = 5

dx = 200
dy = 200

Xpsi = Xfjord  # vertical total
Ypsi = Yfjord  # horizontal total

Xrho = Xpsi+2
Yrho = Ypsi+2
Xvert = Xpsi+3
Yvert = Ypsi+3

# location of fjord mouth
Xm = Xfjord*dx  # m

# fjord width
Yw = Yfjord*dy  # m

# fjord depth
Dm = 200.

# Coriolis Parameter
f0 = 1.e-4
beta = 0.

# vertical grid specs
theta_b = 2.0
theta_s = 8.0
Tcline = 10
N = 20
hsill = 100  # m
xsill = 10000  # m
wsill = 2000  # m

# grid name
grd_name = 'fjord_test'

# ------------ horizontal grid construction ----------------------------
# X direction
xfjord = np.linspace(0, Xm, Xfjord+1)
dxfjord = np.diff(xfjord).mean()

xvert = xfjord
xvert = np.concatenate((np.array([2*xvert[0]-xvert[1]]),
                        xvert,
                        np.array([2*xvert[-1]-xvert[-2]])))

# Y direction
yfjord = np.linspace(-Yw/2., Yw/2., Yfjord+1)
dyfjord = np.diff(yfjord).mean()

yvert = yfjord
yvert = np.concatenate((np.array([2*yvert[0]-yvert[1]]),
                        yvert,
                        np.array([2*yvert[-1]-yvert[-2]])))

# meshgrid it
xxvert, yyvert = np.meshgrid(xvert, yvert)

# generate land mask
msk = np.ones((Yrho, Xrho))
# west, north and south boundary
msk[:, 0] = 0
msk[:, 1] = 0
msk[0, :] = 0
msk[-1, :] = 0

# ------------ write hgrid ---------------------------------------------
hgrd = pyroms.hgrid.CGrid(xxvert, yyvert)
# Coriolis Parameter
hgrd.f = f0 + hgrd.y_rho*beta
hgrd.mask_rho = msk

# ------------ vertical grid construction ------------------------------
h = Dm*np.ones((Yrho, Xrho))
# sill
msks = (hgrd.x_rho >= xsill-wsill) & (hgrd.x_rho <= xsill+wsill)
h[msks] = Dm + (hsill - Dm)*np.exp(-(hgrd.x_rho[msks]-xsill)**2/(0.1*wsill**2))

# ------------ write vgrid ---------------------------------------------
vgrd = pyroms.vgrid.s_coordinate_4(h, theta_b, theta_s, Tcline, N, hraw=h)

# ------------ write grid ----------------------------------------------
grd = pyroms.grid.ROMS_Grid(grd_name, hgrd, vgrd)
pyroms.grid.write_ROMS_grid(grd, filename='/Users/cw686/roms_stuff/grid/fjord_grd_test.nc')

