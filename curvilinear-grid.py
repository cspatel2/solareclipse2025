#%%
import matplotlib.pyplot as plt
import numpy as np

from matplotlib.projections import PolarAxes
from matplotlib.transforms import Affine2D
from mpl_toolkits.axisartist import Axes, HostAxes, angle_helper
from mpl_toolkits.axisartist.grid_helper_curvelinear import \
    GridHelperCurveLinear

#%%
def curvelinear_test2(fig):
    """
    Polar projection, but in a rectangular box.
    """

    # PolarAxes.PolarTransform takes radian. However, we want our coordinate
    # system in degree
    tr = Affine2D().scale(np.pi/180, np.pi/180) + PolarAxes.PolarTransform(apply_theta_transforms=False)
    # Polar projection, which involves cycle, and also has limits in
    # its coordinates, needs a special method to find the extremes
    # (min, max of the coordinate within the view).
    extreme_finder = angle_helper.ExtremeFinderCycle(
        nx=20, ny=20,  # Number of sampling points in each direction.
        lon_cycle=360, lat_cycle=None,
        lon_minmax=(0,360), lat_minmax=(0,90 ),
    )
    # Find grid values appropriate for the coordinate (degree, minute, second).
    grid_locator1 = angle_helper.LocatorDMS(5)
    # Use an appropriate formatter.  Note that the acceptable Locator and
    # Formatter classes are a bit different than that of Matplotlib, which
    # cannot directly be used here (this may be possible in the future).
    tick_formatter1 = angle_helper.FormatterDMS()

    grid_helper = GridHelperCurveLinear(
        tr, 
        extreme_finder=extreme_finder,
        grid_locator1=grid_locator1, tick_formatter1=tick_formatter1)
    ax1 = fig.add_subplot(
        1, 2, 2, axes_class=HostAxes, grid_helper=grid_helper)

    # make ticklabels of right and top axis visible.
    ax1.axis["right"].major_ticklabels.set_visible(True)
    ax1.axis["top"].major_ticklabels.set_visible(True)
    # let right axis shows ticklabels for 1st coordinate (angle)
    ax1.axis["right"].get_helper().nth_coord_ticks = 0
    # let bottom axis shows ticklabels for 2nd coordinate (radius)
    ax1.axis["bottom"].get_helper().nth_coord_ticks = 1

    ax1.set_aspect(1)
    ax1.set_xlim(-1, 1)
    ax1.set_ylim(-1, 1)

    ax1.grid(True)

    # # A parasite Axes with given transform
    # ax2 = ax1.get_aux_axes(tr)
    # # note that ax2.transData == tr + ax1.transData
    # # Anything you draw in ax2 will match the ticks and grids of ax1.
    # ax2.plot(np.linspace(0, 30, 51), np.linspace(10, 10, 51), linewidth=2)

    # ax2.pcolor(np.linspace(0, 90, 4), np.linspace(0, 10, 4),
    #            np.arange(9).reshape((3, 3)))
    # ax2.contour(np.linspace(0, 90, 4), np.linspace(0, 10, 4),
    #             np.arange(16).reshape((4, 4)), colors="k")


if __name__ == "__main__":
    fig = plt.figure(figsize=(7, 4))
    curvelinear_test2(fig)

    plt.show()






# %%
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
from skyfield.api import Star, load
from skyfield.data import hipparcos, stellarium
from skyfield.projections import build_stereographic_projection
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs import WCS
from astropy.visualization.wcsaxes import WCSAxes
import numpy as np

from matplotlib.projections import PolarAxes
from matplotlib.transforms import Affine2D
from mpl_toolkits.axisartist import Axes, HostAxes, angle_helper
from mpl_toolkits.axisartist.grid_helper_curvelinear import \
    GridHelperCurveLinear
# %%
AZ = 230
ALT = 75
COORD = SkyCoord(f'{AZ} {ALT}', unit = (u.deg, u.deg), frame='altaz')
fig = plt.figure(figsize=[6,6])
# wcs = WCS(naxis = 2)
# wcs.wcs.crpix = [1,1]
# wcs.wcs.crpix = [1, 1]
# wcs.wcs.cdelt = np.array([-360 / np.pi, 360/np.pi])
# wcs.wcs.crval = [COORD.az.deg, COORD.alt.deg]
# # wcs.wcs.ctype = ["RA---STG", "DEC--STG"]
# print(f'Alt, Az = {COORD.alt.deg},{COORD.az.deg}')

tr = Affine2D().scale(np.pi/180, 90) + PolarAxes.PolarTransform(apply_theta_transforms=False)
grid_locator1 = angle_helper.LocatorDMS(12)

grid_helper = GridHelperCurveLinear(
    tr, grid_locator1=grid_locator1)

ax = fig.add_subplot(
    1, 1, 1, axes_class=HostAxes, grid_helper=grid_helper)
ax.axis["right"].major_ticklabels.set_visible(True)
ax.axis["top"].major_ticklabels.set_visible(True)
ax.axis["right"].get_helper().nth_coord_ticks = 0
ax.axis["bottom"].get_helper().nth_coord_ticks = 1

ax.set_aspect(1)
ax.grid(True, zorder=0)
ax.grid(True, color='black', linestyle='dotted')
# <- End point code implementetion from docs




# %%
