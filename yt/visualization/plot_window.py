"""
A plotting mechanism based on the idea of a "window" into the data.

Author: J. S. Oishi <jsoishi@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2010-2011 J. S. Oishi.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
import base64
import matplotlib.pyplot
from functools import wraps

import numpy as na
from .image_writer import \
    write_image, apply_colormap
from .fixed_resolution import \
    FixedResolutionBuffer
from .plot_modifications import get_smallest_appropriate_unit
from .tick_locators import LogLocator, LinearLocator

from yt.funcs import *
from yt.utilities.amr_utils import write_png_to_string

def invalidate_data(f):
    @wraps(f)
    def newfunc(*args, **kwargs):
        f(*args, **kwargs)
        args[0]._data_valid = False
        args[0]._plot_valid = False
        args[0]._recreate_frb()
        if args[0]._initfinished:
            args[0]._setup_plots()
    return newfunc

def invalidate_plot(f):
    @wraps(f)
    def newfunc(*args, **kwargs):
        args[0]._plot_valid = False
        args[0]._setup_plots()
        return f(*args, **kwargs)
    return newfunc

class PlotWindow(object):
    def __init__(self, data_source, bounds, buff_size=(800,800), antialias = True, periodic = True):
        r"""
        PlotWindow(data_source, bounds, buff_size=(800,800), antialias = True)
        
        A ploting mechanism based around the concept of a window into a
        data source. It can have arbitrary fields, each of which will be
        centered on the same viewpoint, but will have individual zlimits. 
        
        The data and plot are updated separately, and each can be
        invalidated as the object is modified.
        
        Data is handled by a FixedResolutionBuffer object.

        Parameters
        ----------
        data_source : :class:`yt.data_objects.data_containers.AMRProjBase` or :class:`yt.data_objects.data_containers.AMRSliceBase`
            This is the source to be pixelized, which can be a projection or a
            slice.  (For cutting planes, see
            `yt.visualization.fixed_resolution.ObliqueFixedResolutionBuffer`.)
        bounds : sequence of floats
            Bounds are the min and max in the image plane that we want our
            image to cover.  It's in the order of (xmin, xmax, ymin, ymax),
            where the coordinates are all in the appropriate code units.
        buff_size : sequence of ints
            The size of the image to generate.
        antialias : boolean
            This can be true or false.  It determines whether or not sub-pixel
            rendering is used during data deposition.

        """
        self._initfinished = False
        self.center = None
        self.plots = {}
        self._periodic = periodic
        self.data_source = data_source
        self.buff_size = buff_size
        self.antialias = True
        self.set_window(bounds) # this automatically updates the data and plot
        if self.data_source.center is not None:
            center = [self.data_source.center[i] for i in range(len(self.data_source.center)) if i != self.data_source.axis]
            self.set_center(center)
        self._initfinished = True

    def __getitem__(self, item):
        return self.plots[item]

    def _recreate_frb(self):
        try:
            bounds = self.bounds
            self._frb = FixedResolutionBuffer(self.data_source, 
                                              bounds, self.buff_size, 
                                              self.antialias, periodic=self._periodic)
        except:
            raise RuntimeError("Failed to repixelize.")
        self._frb._get_data_source_fields()
        self._data_valid = True
        
    def _setup_plots(self):
        pass

    @property
    def fields(self):
        return self._frb.data.keys()

    @property
    def width(self):
        Wx = self.xlim[1] - self.xlim[0]
        Wy = self.ylim[1] - self.ylim[0]
        return (Wx, Wy)

    @property
    def bounds(self):
        return self.xlim+self.ylim

    @invalidate_data
    def zoom(self, factor):
        r"""This zooms the window by *factor*.

        Parameters
        ----------
        factor : float
            multiplier for the current width

        """
        Wx, Wy = self.width
        centerx = self.xlim[0] + Wx*0.5
        centery = self.ylim[0] + Wy*0.5
        nWx, nWy = Wx/factor, Wy/factor
        self.xlim = (centerx - nWx*0.5, centerx + nWx*0.5)
        self.ylim = (centery - nWy*0.5, centery + nWy*0.5)

    @invalidate_data
    def pan(self, deltas):
        r"""Pan the image by specifying absolute code unit coordinate deltas.
        
        Parameters
        ----------
        deltas : sequence of floats
            (delta_x, delta_y) in *absolute* code unit coordinates

        """
        self.xlim = (self.xlim[0] + deltas[0], self.xlim[1] + deltas[0])
        self.ylim = (self.ylim[0] + deltas[1], self.ylim[1] + deltas[1])

    @invalidate_data
    def pan_rel(self, deltas):
        r"""Pan the image by specifying relative deltas, to the FOV.
        
        Parameters
        ----------
        deltas : sequence of floats
            (delta_x, delta_y) in *relative* code unit coordinates

        """
        Wx, Wy = self.width
        self.xlim = (self.xlim[0] + Wx*deltas[0], self.xlim[1] + Wx*deltas[0])
        self.ylim = (self.ylim[0] + Wy*deltas[1], self.ylim[1] + Wy*deltas[1])

    @invalidate_data
    def set_field(self):
        pass

    @invalidate_data
    def set_window(self, bounds):
        if self.center is not None:
            dx = bounds[1] - bounds[0]
            dy = bounds[3] - bounds[2]
            self.xlim = (self.center[0] - dx/2., self.center[0] + dx/2.)
            self.ylim = (self.center[1] - dy/2., self.center[1] + dy/2.)
            mylog.info("xlim = %f %f" %self.xlim)
            mylog.info("ylim = %f %f" %self.ylim)
        else:
            self.xlim = bounds[0:2]
            self.ylim = bounds[2:]
        
    @invalidate_data
    def set_width(self, new_width):
        """set the width of the plot window

        parameters
        ----------
        new_width : float
            the width of the image in code units.

        """
        Wx, Wy = self.width
        centerx = self.xlim[0] + Wx*0.5
        centery = self.ylim[0] + Wy*0.5
        self.xlim = (centerx - new_width/2.,
                     centerx + new_width/2.)
        self.ylim = (centery - new_width/2.,
                     centery + new_width/2.)

    @invalidate_data
    def set_center(self, new_center):
        if new_center is None:
            self.center = None
        else:
            self.center = new_center
        self.set_window(self.bounds)

    @property
    def width(self):
        Wx = self.xlim[1] - self.xlim[0]
        Wy = self.ylim[1] - self.ylim[0]
        return (Wx, Wy)

    # @property
    # def window(self):
    #     return self.xlim + self.ylim

    @invalidate_data
    def set_antialias(self,aa):
        self.antialias = aa

class PWViewer(PlotWindow):
    """A viewer for PlotWindows.

    """
    def __init__(self, *args,**kwargs):
        setup = kwargs.pop("setup", True)
        PlotWindow.__init__(self, *args,**kwargs)
        self._field_transform = {}
        for field in self._frb.data.keys():
            if self._frb.pf.field_info[field].take_log:
                self._field_transform[field] = na.log
            else:
                self._field_transform[field] = lambda x: x

        if setup: self._setup_plots()

    @invalidate_plot
    def set_log(self,field,log):
        """set a field to log or linear.
        
        Parameters
        ----------
        field : string
            the field to set a transform
        log : boolean
            Log on/off.

        """
        if log:
            self._field_transform[field] = na.log
        else:
            self._field_transform[field] = lambda x: x

    def set_transform(self, field, func):
        self._field_transform[field] = func

    @invalidate_plot
    def set_cmap(self):
        pass

    @invalidate_plot
    def set_zlim(self):
        pass

class PWViewerMPL(PWViewer):
    """Viewer using matplotlib as a backend via the YtWindowPlot. 

    """
    def _setup_plots(self):
        for f in self.fields:
            self.plots[f] = YtWindowPlot(self._frb[f])
        self._plot_valid = True

    def save(self,name):
        for k,v in self.plots.iteritems():
            n = "%s_%s" % (name, k)
            v.save(n)

class PWViewerRaw(PWViewer):
    """A PlotWindow viewer that writes raw pngs (no MPL, no axes).

    """
    def _setup_plots(self):
        self.save('')
        self._plot_valid = True

    def save(self,name):
        for field in self._frb.data.keys():
            nm = "%s_%s.png" % (name,field)
            print "writing %s" % nm
            write_image(self._frb[field],nm)

_metadata_template = """
%(pf)s<br>
<br>
Field of View:  %(x_width)0.3f %(unit)s<br>
Minimum Value:  %(mi)0.3e %(units)s<br>
Maximum Value:  %(ma)0.3e %(units)s
"""

class PWViewerExtJS(PWViewer):
    """A viewer for the web interface.

    """
    _ext_widget_id = None
    _current_field = None
    _widget_name = "plot_window"
    cmap = 'algae'

    def _setup_plots(self):
        from yt.gui.reason.bottle_mods import PayloadHandler
        ph = PayloadHandler()
        if self._current_field is not None \
           and self._ext_widget_id is not None:
            fields = [self._current_field]
            addl_keys = {'type': 'widget_payload',
                         'widget_id': self._ext_widget_id}
        else:
            fields = self._frb.data.keys()
            addl_keys = {}
        min_zoom = 200*self._frb.pf.h.get_smallest_dx() * self._frb.pf['unitary']
        for field in fields:
            to_plot = apply_colormap(self._frb[field], func = self._field_transform[field])
            pngs = write_png_to_string(to_plot)
            img_data = base64.b64encode(pngs)
            # We scale the width between 200*min_dx and 1.0
            x_width = self.xlim[1] - self.xlim[0]
            zoom_fac = na.log10(x_width*self._frb.pf['unitary'])/na.log10(min_zoom)
            zoom_fac = 100.0*max(0.0, zoom_fac)
            ticks = self.get_ticks(self._frb[field].min(),
                                   self._frb[field].max(), 
                                   take_log = self._frb.pf.field_info[field].take_log)
            payload = {'type':'png_string',
                       'image_data':img_data,
                       'metadata_string': self.get_metadata(field),
                       'zoom': zoom_fac,
                       'ticks': ticks}
            payload.update(addl_keys)
            ph.add_payload(payload)

    def get_ticks(self, mi, ma, height = 400, take_log = False):
        # This will eventually change to work with non-logged fields
        ticks = []
        if take_log:
            ll = LogLocator() 
            tick_locs = ll(mi, ma)
            mi = na.log10(mi)
            ma = na.log10(ma)
            for v1,v2 in zip(tick_locs, na.log10(tick_locs)):
                if v2 < mi or v2 > ma: continue
                p = height - height * (v2 - mi)/(ma - mi)
                ticks.append((p,v1,v2))
                #print v1, v2, mi, ma, height, p
        else:
            ll = LinearLocator()
            tick_locs = ll(mi, ma)
            for v in tick_locs:
                p = height - height * (v - mi)/(ma-mi)
                ticks.append((p,v,"%0.3e" % (v)))

        return ticks

    def _get_cbar_image(self, height = 400, width = 40):
        # Right now there's just the single 'cmap', but that will eventually
        # change.  I think?
        vals = na.mgrid[1:0:height * 1j] * na.ones(width)[:,None]
        vals = vals.transpose()
        to_plot = apply_colormap(vals)
        pngs = write_png_to_string(to_plot)
        img_data = base64.b64encode(pngs)
        return img_data

    # This calls an invalidation routine from within
    def scroll_zoom(self, value):
        # We accept value from 0..100, and assume it has been set from the
        # scroll bar.  In that case, we undo the logic for calcualting
        # 'zoom_fac' from above.
        min_val = 200*self._frb.pf.h.get_smallest_dx()
        unit = self._frb.pf['unitary']
        width = (min_val**(value/100.0))/unit
        self.set_width(width)

    def get_metadata(self, field):
        fval = self._frb[field]
        mi = fval.min()
        ma = fval.max()
        x_width = self.xlim[1] - self.xlim[0]
        y_width = self.ylim[1] - self.ylim[0]
        unit = get_smallest_appropriate_unit(x_width, self._frb.pf)
        units = self.get_field_units(field)
        md = _metadata_template % dict(
                pf = self._frb.pf,
                x_width = x_width*self._frb.pf[unit],
                y_width = y_width*self._frb.pf[unit],
                unit = unit, units = units, mi = mi, ma = ma)
        return md

    def image_recenter(self, img_x, img_y, img_size_x, img_size_y):
        dx = (self.xlim[1] - self.xlim[0]) / img_size_x
        dy = (self.ylim[1] - self.ylim[0]) / img_size_y
        new_x = img_x * dx + self.xlim[0]
        new_y = img_y * dy + self.ylim[0]
        print img_x, img_y, dx, dy, new_x, new_y
        self.set_center((new_x, new_y))

    @invalidate_data
    def set_current_field(self, field):
        self._current_field = field
        self._frb[field]
        if self._frb.pf.field_info[field].take_log:
            self._field_transform[field] = na.log
        else:
            self._field_transform[field] = lambda x: x

    def get_field_units(self, field, strip_mathml = True):
        ds = self._frb.data_source
        pf = self._frb.pf
        if ds._type_name == "slice":
            units = pf.field_info[field].get_units()
        elif ds._type_name == "proj":
            units = pf.field_info[field].get_projected_units()
        else:
            units = ""
        if strip_mathml:
            units = units.replace(r"\rm{", "").replace("}","")
        return units


class YtPlot(object):
    """A base class for all yt plots. It should abstract the actual
    plotting engine completely, allowing plotting even without matplotlib. 

    YtPlot and the classes that derive from it are *by design* limited
    and designed for rapid, production quality plot production, rather
    than full generality. If you require more customization of the end
    result, these objects are designed to return to you the basic data
    so you the user can insert them into a matplotlib figure on your
    own outside of the YtPlot class.

    """
    axis_names = {}
    datalabel = None
    figure = None
    def __init__(self, field, size=(10,8)):
        self.__setup_from_field(field)
        self._plot_valid = True
        self.figure = matplotlib.pyplot.figure(figsize=size)
        self.axes = self.figure.add_subplot(1,1,1)

    def __setup_from_field(self, field):
        #self.set_log_field(self.pf.field_info[field].take_log)
        self.axis_names["X"] = None
        self.axis_names["Y"] = None
        self.axis_names["Z"] = field

    def save(self,name):
        print "saving plot %s" % name
        self.figure.savefig('%s.png' % name)

class Yt2DPlot(YtPlot):
    zmin = None
    zmax = None
    cmap = 'algae'
    zlabel = None

    # def __init__(self, data):
    #     pass

    @invalidate_plot
    def set_zlim(self, zmin, zmax):
        self.zmin = zmin
        self.zmax = zmax

    @invalidate_plot
    def set_cmap(self,cmap):
        self.cmap = cmap

class YtWindowPlot(Yt2DPlot):
    def __init__(self, data, size=(10,8)):
        YtPlot.__init__(self, data, size)
        self.__init_image(data)

    def __init_image(self, data):
        self.image = self.axes.imshow(data,cmap=self.cmap)

class YtProfilePlot(Yt2DPlot):
    def __init__(self):
        pass