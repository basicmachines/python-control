import numpy as np
from numpy import cos, sin, sqrt, linspace, pi, exp
import matplotlib.pyplot as plt
from mpl_toolkits.axisartist import SubplotHost, axislines
from mpl_toolkits.axisartist.grid_helper_curvelinear import GridHelperCurveLinear
import mpl_toolkits.axisartist.angle_helper as angle_helper
from matplotlib.projections import PolarAxes
from matplotlib.transforms import Affine2D


class FormatterDMS(object):
    """Transforms angle ticks to damping ratios"""

    def __call__(self, direction, factor, values):
        angles_deg = values/factor
        damping_ratios = np.cos((180-angles_deg)*np.pi/180)
        ret = ["%.2f" % val for val in damping_ratios]
        return ret


class ModifiedExtremeFinderCycle(angle_helper.ExtremeFinderCycle):
    """Changed to allow only left hand-side polar grid"""

    def __call__(self, transform_xy, x1, y1, x2, y2):
        x_, y_ = np.linspace(x1, x2, self.nx), np.linspace(y1, y2, self.ny)
        x, y = np.meshgrid(x_, y_)
        lon, lat = transform_xy(np.ravel(x), np.ravel(y))

        with np.errstate(invalid='ignore'):
            if self.lon_cycle is not None:
                lon0 = np.nanmin(lon)
                # Changed from 180 to 360 to be able to span only 90-270 (left hand side)
                lon -= 360. * ((lon - lon0) > 360.)
            if self.lat_cycle is not None:
                lat0 = np.nanmin(lat)
                # Changed from 180 to 360 to be able to span only 90-270 (left hand side)
                lat -= 360. * ((lat - lat0) > 360.)

        lon_min, lon_max = np.nanmin(lon), np.nanmax(lon)
        lat_min, lat_max = np.nanmin(lat), np.nanmax(lat)

        lon_min, lon_max, lat_min, lat_max = \
            self._adjust_extremes(lon_min, lon_max, lat_min, lat_max)

        return lon_min, lon_max, lat_min, lat_max


def sgrid(ax=None):
    """
    Adds an s-plane grid of constant damping factors and natural
    frequencies to a plot. If ax is not specified, the current
    figure axis is used. However, note that the returned axis
    may not be the same instance as the axis passed. The original
    axis may have to be deleted to create the s-plane axis.
    Always use the new axis instance returned by this method if
    you need to reference the axis created.

    Parameters
    ----------
    ax : matplotlib axis
        If not passed, uses gca() to get current axis.

    Returns
    -------
    ax : matplotlib axis

    Example
    -------
    >>> H = tf([2, 5, 1], [1, 2, 3])
    >>> r, k = rlocus(H)
    >>> ax = sgrid()
    >>> plt.show()
    """
    # Based on these matplotlib demos:
    # https://matplotlib.org/gallery/axisartist/demo_curvelinear_grid.html
    # https://matplotlib.org/gallery/axisartist/demo_floating_axis.html

    # PolarAxes.PolarTransform takes radian. However, we want our
    # coordinate system in degrees
    tr = Affine2D().scale(np.pi/180., 1.) + PolarAxes.PolarTransform()
    # polar projection, which involves cycle, and also has limits in
    # its coordinates, needs a special method to find the extremes
    # (min, max of the coordinate within the view).

    # 20, 20 : number of sampling points along x, y direction
    sampling_points = 20
    extreme_finder = ModifiedExtremeFinderCycle(sampling_points,
                                                sampling_points,
                                                lon_cycle=360,
                                                lat_cycle=None,
                                                lon_minmax=(90, 270),
                                                lat_minmax=(0, np.inf))
    grid_locator1 = angle_helper.LocatorDMS(15)
    tick_formatter1 = FormatterDMS()
    grid_helper = GridHelperCurveLinear(tr,
                                        extreme_finder=extreme_finder,
                                        grid_locator1=grid_locator1,
                                        tick_formatter1=tick_formatter1
                                        )

    if ax is None:
        # Get the current axis or create a new figure with one
        ax = plt.gca()
    fig = ax.figure

    # If the axis is not from mpl_toolkits.axisartist, replace it
    # with a new axis
    if not isinstance(ax, axislines.Axes):
        subargs = ax.get_geometry()
        fig.delaxes(ax)
        # TODO: The caller of sgrid() still references ax which is no
        #       longer part of the figure. Is there a more sophisticated
        #       way to replace the axes reference?
        ax2 = SubplotHost(fig, *subargs, grid_helper=grid_helper)
        fig.add_subplot(ax2)
    else:
        ax2 = ax

    # make tick labels of right invisible, and top axis visible.
    visible = True
    ax2.axis[:].major_ticklabels.set_visible(visible)
    ax2.axis[:].major_ticks.set_visible(False)
    ax2.axis[:].invert_ticklabel_direction()

    ax2.axis["wnxneg"] = axis = ax2.new_floating_axis(0, 180)
    axis.set_ticklabel_direction("-")
    axis.label.set_visible(False)
    ax2.axis["wnxpos"] = axis = ax2.new_floating_axis(0, 0)
    axis.label.set_visible(False)
    ax2.axis["wnypos"] = axis = ax2.new_floating_axis(0, 90)
    axis.label.set_visible(False)
    axis.set_axis_direction("left")
    ax2.axis["wnyneg"] = axis = ax2.new_floating_axis(0, 270)
    axis.label.set_visible(False)
    axis.set_axis_direction("left")
    axis.invert_ticklabel_direction()
    axis.set_ticklabel_direction("-")

    # let left axis shows ticklabels for 1st coordinate (angle)
    ax2.axis["left"].get_helper().nth_coord_ticks = 0
    ax2.axis["right"].get_helper().nth_coord_ticks = 0
    ax2.axis["left"].get_helper().nth_coord_ticks = 0
    ax2.axis["bottom"].get_helper().nth_coord_ticks = 0

    # TODO: Not sure how to add this axis
    #fig.add_subplot(ax)

    # RECTANGULAR X Y AXES WITH SCALE
    #par2 = ax.twiny()
    #par2.axis["top"].toggle(all=False)
    #par2.axis["right"].toggle(all=False)
    #new_fixed_axis = par2.get_grid_helper().new_fixed_axis
    #par2.axis["left"] = new_fixed_axis(loc="left",
    #                                   axes=par2,
    #                                   offset=(0, 0))
    #par2.axis["bottom"] = new_fixed_axis(loc="bottom",
    #                                     axes=par2,
    #                                     offset=(0, 0))
    # FINISH RECTANGULAR

    ax2.grid(True, zorder=0, linestyle='dotted')
    _final_setup(ax2)
    return ax2


def _final_setup(ax):
    ax.set_xlabel('Real')
    ax.set_ylabel('Imaginary')
    ax.axhline(y=0, color='black', lw=1)
    ax.axvline(x=0, color='black', lw=1)
    plt.axis('equal')


def nogrid(ax=None):
    if ax is None:
        # Get the current axis or create a new figure with one
        ax = plt.gca()

    _final_setup(ax)
    return ax


def zgrid(ax=None, zetas=None, wns=None):
    """Adds a z-plane grid of constant damping factors and natural
    frequencies to a plot. If ax is not specified, the current
    figure axis is used.

    Parameters
    ----------
    ax : matplotlib axis
        If not passed, uses gca() to get current axis.
    zetas : list or 1-d array
        Damping factors.
    wns : list or 1-d array
        Normalized natural frequencies.

    Returns
    -------
    ax : matplotlib axis

    Example
    -------
    >>> H = tf([2, -3.4, 1.5], [1, -1.6, 0.8], -1)
    >>> r, k = rlocus(H)
    >>> zgrid()
    >>> plt.show()
    """

    if ax is None:
        # Get the current axis or create a new figure with one
        ax = plt.gca()

    # Constant damping lines
    if zetas is None:
        zetas = linspace(0, 0.9, 10)
    for zeta in zetas:
        # Calculate in polar coordinates
        factor = zeta/sqrt(1-zeta**2)
        x = linspace(0, sqrt(1-zeta**2), 200)
        ang = pi*x
        mag = exp(-pi*factor*x)
        # Draw upper part in rectangular coordinates
        xret = mag*cos(ang)
        yret = mag*sin(ang)
        ax.plot(xret, yret, 'k:', lw=1)
        # Draw lower part in rectangular coordinates
        xret = mag*cos(-ang)
        yret = mag*sin(-ang)
        ax.plot(xret, yret, 'k:', lw=1)
        # Annotation
        an_i = int(len(xret)/2.5)
        an_x = xret[an_i]
        an_y = yret[an_i]
        ax.annotate(str(round(zeta, 2)), xy=(an_x, an_y),
                    xytext=(an_x, an_y), size=7)

    # Constant natural frequency lines
    if wns is None:
        wns = linspace(0, 1, 10)
    for a in wns:
        # Calculate in polar coordinates
        x = linspace(-pi/2, pi/2, 200)
        ang = pi*a*sin(x)
        mag = exp(-pi*a*cos(x))
        # Draw in rectangular coordinates
        xret = mag*cos(ang)
        yret = mag*sin(ang)
        ax.plot(xret, yret, 'k:', lw=1)
        # Annotation
        an_i = -1
        an_x = xret[an_i]
        an_y = yret[an_i]
        num = '{:1.1f}'.format(a)
        ax.annotate("$\\frac{" + num + "\\pi}{T}$", xy=(an_x, an_y),
                    xytext=(an_x, an_y), size=9)

    _final_setup(ax)
    return ax
