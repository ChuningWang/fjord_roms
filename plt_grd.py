import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

from cmocean import cm
import pyroms

grd = pyroms.grid.get_ROMS_grid('fjord_test')
x = grd.hgrid.x_vert
y = grd.hgrid.y_vert
mask = grd.hgrid.mask_rho
mask[1:-1, 1] = 2
h = grd.vgrid.h

cmap = colors.ListedColormap(['k', 'w', 'b'])

fig, ax = plt.subplots()
ax.scatter(-500, 0,
           marker='s', facecolor='k', edgecolor='k',
           alpha=0.5, label='Land')
# ax.scatter(-500, 0,
#            marker='s', facecolor='w', edgecolor='k',
#            alpha=0.5, label='Ocean')
ax.scatter(-500, 0,
           marker='s', facecolor='b', edgecolor='k',
           alpha=0.5, label='Glacier')
ax.set_aspect('equal')
ax.set_xlim(-200, 10200)
ax.set_ylim(-700, 700)
ax.set_xlabel('X [m]')
ax.set_ylabel('Y [m]')
ax.plot(x.T, y.T, 'k')
ax.plot(x, y, 'k')
ax.pcolor(x, y, mask, alpha=0.5, cmap=cmap)
ax.legend(loc=9, ncol=3, bbox_to_anchor=(0.25, 1.6))
pc = ax.pcolor(x[1:-1, 2:], y[1:-1, 2:], h[1:-1, 2:], vmin=50, vmax=200)
cbar_ax = fig.add_axes([0.52, 0.58, 0.35, 0.02])
cb = fig.colorbar(pc, cax=cbar_ax,
        ticks=np.linspace(50, 200, 7), orientation='horizontal')
cb.ax.tick_params(labelsize=5)
cb.set_label('Depth [m]', fontsize=7)
cbar_ax.xaxis.tick_top()
cbar_ax.xaxis.set_label_position('top')
fig.savefig('grid.png', dpi=300)
plt.close()
