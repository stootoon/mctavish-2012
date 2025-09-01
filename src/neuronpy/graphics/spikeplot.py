# -*- coding: utf-8 -*-
"""Plot spike data in a raster plot.

Plot spike times as a raster plot where the horizontal axis is time and the
rows of the vertical axis are spike trains of individual neurons. This 
sets and holds many variables for plotting such that the user
can simply set the parameters once and then redraw the spikes or scale and
zoom to navigate spike data.

EXAMPLE:

This simple example generates some random "spikes" and displays them.

::

    import numpy
    from neuronpy.graphics import spikeplot

    spikes = []
    num_cells = 10
    num_spikes_per_cell = 20
    frequency = 20

    # Make the spike data. Use a simple Poisson-like spike generator 
    # (just for illustrative purposes here. Better spike generators should 
    # be used in simulations).
    for i in range(num_cells):
        isi=numpy.random.poisson(frequency, num_spikes_per_cell)
        spikes.append(numpy.cumsum(isi))
        
    # spikes is now a list of lists where each cell has a list of spike
    # times. Now, let's plot these spikes.
    sp = spikeplot.SpikePlot(savefig=True)
    sp.plot_spikes(spikes)

AUTHORS:

- THOMAS MCTAVISH (2010-03-01): initial version, 0.1
- THOMAS MCTAVISH (2012-04-10): Few changes to get spiketimehistograms and
    spikesums to coexist better.
"""
# While this software is under the permissive MIT License, 
# (http://www.opensource.org/licenses/mit-license.php)
# We ask that you cite the neuronpy package (or tools used in this package)
# in any publications and contact the author with your referenced publication.
#
# Format:
# McTavish, T.S. NeuronPy library, version 0.1.5, http://bitbucket.org/tommctavish/neuronpy
#
# Copyright (c) 2012 Thomas S. McTavish
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import sys

import numpy
from matplotlib import pyplot
from matplotlib import lines
from matplotlib.colors import colorConverter
from matplotlib.backends.backend_pdf import PdfPages

from neuronpy.graphics import spikeparams
from neuronpy.graphics import spiketimehistogram
from neuronpy.graphics import spikesum
from neuronpy.util import listutil
from neuronpy.util import spiketrain
from neuronpy.util import vartest
from neuronpy.util.decorators import set_overrides_and_defaults
    
__version__ = "0.1.1"
            
class SpikePlot(object):
    """
    Manager for the spike raster 
    `axes <http://matplotlib.sourceforge.net/api/axes_api.html>`_ as well 
    as the axes manager for 
    :class:`~neuronpy.graphics.spiketimehistogram.SpikeTimeHistogram`, and
    :class:`~neuronpy.graphics.spikesum.SpikeSum` subplots. Additionally, the
    line formatting parameters in 
    :class:`~neuronpy.graphics.spikeparams.SpikeParams`
    funnel through this object. In fact, it may be rare to access those
    classes directly and instead use the methods in this class
    
    Plots are drawn iteratively in that the user specifies the drawing 
    parameters for a set of spikes *before* telling the ``SpikePlot`` object to
    plot any spikes. Therefore, when plotting multiple sets of spikes in the
    same plot, the paradigm is::
    
        sp = SpikePlot()
        sp.plot_spikes(some_spikes)
        sp.set_<some_line_format_method>()
        sp.set_<some_other_formatting_options>()
        sp.plot_spikes(some_other_spikes)
        
    Indeed, the method to be most familiar with is :func:`plot_spikes`.
    
    If the spike data is not too large to fit into RAM and you want to 
    navigate it, it is generally a good idea to specify the complete data 
    to plot and then simply modify the view of the axes to zoom or translate 
    the plot. This is detailed in :func:`update_xlim`.
    
    Parameters can be set during instanciation of a ``SpikePlot`` object by
    passing in various ``<parameter>=<value>`` pairs or through accessor 
    methods, ``get_<parameter>()`` and ``set_<parameter>()``. As 
    ``<parameter>=<value>`` pairs, they use the ``kwargs`` (short for 
    *keyword arguments*) parameter on instanciation. Below is a summary of 
    keyword parameters, which are more detailed in the corresponding 
    ``set_<parameter>()`` method. Some lesser-used parameters have accessor
    methods, but are not available as keyword parameters.

    Summary of SpikePlot keyword parameters:
    
    +------------------+---------------------------------+-------------------+
    |     Parameter    |            Description          |      Default      |
    +==================+=================================+===================+
    | fig_             | The Figure_ handle.             | ``None``          |
    +------------------+---------------------------------+-------------------+
    | fig_name_        | Name of figure, if saved as     |``'spikeplot.png'``|
    |                  | a file.                         |                   |
    +------------------+---------------------------------+-------------------+
    | figsize_         | Size of the figure, in inches.  | ``(8, 6)``        |
    +------------------+---------------------------------+-------------------+
    | marker_          | Tick marks style.               | ``'|'``           |
    +------------------+---------------------------------+-------------------+
    | markercolor_     | Line and marker color.          | ``'black'``       |
    +------------------+---------------------------------+-------------------+
    | markeredgewidth_ | Line and edge width of tick     | ``0.75``          |
    |                  | marks.                          |                   |
    +------------------+---------------------------------+-------------------+
    | markerscale_     | Tick mark row height.           | ``1.0``           |
    +------------------+---------------------------------+-------------------+
    | savefig_         | Flag to save a figure to a      | ``False``         |
    |                  | file when drawing.              |                   |
    +------------------+---------------------------------+-------------------+
    | sth_ratio_       | |SpikeTimeHistogram| axes       | ``0``             |
    |                  | height.                         |                   |
    +------------------+---------------------------------+-------------------+ 
    | sum_ratio_       | |SpikeSum| Histogram axes       | ``0``             |
    |                  | width.                          |                   |
    +------------------+---------------------------------+-------------------+ 
    
    .. _fig: #neuronpy.graphics.spikeplot.SpikePlot.set_fig
    .. _Figure: http://matplotlib.sourceforge.net/api/figure_api.html
    .. _fig_name: #neuronpy.graphics.spikeplot.SpikePlot.set_fig_name
    .. _figsize: #neuronpy.graphics.spikeplot.SpikePlot.set_figsize
    .. _marker: #neuronpy.graphics.spikeplot.SpikePlot.set_marker
    .. _markercolor: #neuronpy.graphics.spikeplot.SpikePlot.set_markercolor
    .. _markeredgewidth: \
    #neuronpy.graphics.spikeplot.SpikePlot.set_markeredgewidth
    .. _markerscale: #neuronpy.graphics.spikeplot.SpikePlot.set_markerscale
    .. _savefig: #neuronpy.graphics.spikeplot.SpikePlot.set_savefig
    .. _sth_ratio: #neuronpy.graphics.spikeplot.SpikePlot.set_sth_ratio
    .. _sum_ratio: #neuronpy.graphics.spikeplot.SpikePlot.set_sum_ratio
    .. |SpikeSum| replace:: 
        :class:`~neuronpy.graphics.spikesum.SpikeSum`
    .. |SpikeTimeHistogram| replace:: 
        :class:`~neuronpy.graphics.spiketimehistogram.SpikeTimeHistogram`
    """
    def __init__(self, **kwargs):
        """
        Initialize a SpikePlot instance.
        """
        # Make a scratch variable for line formatting tests.
        self._line2d = lines.Line2D([1], [1])
        
        self._raster_axes = None
        self._sum = None
        self._sth = None

        # Make a dictionary of parameters and their default values.
        defaults = { 'fig':None,
                'fig_name':'spikeplot.png',
                'figsize':None,
                'linestyle':'None',
                'linewidth':0.75,
                'marker':'|',
                'markercolor':'black',
                'markeredgewidth':0.75,
                'markerscale':1.0,
                'savefig':False,
                'sth_ratio':0.,
                'sum_ratio':0.}
        
        set_overrides_and_defaults(self, defaults, kwargs)
                
        # Reset the figsize variable, if a figure was passed in
        if self._fig is not None:
            self.set_figsize(self._fig.get_size_inches())
            
        # Set other member variables not in the keywords
        self._data_xlim = [sys.maxsize, 0]
        self._axes_lim = [0, 0]
        self._spike_params = dict()
        self._sth = spiketimehistogram.SpikeTimeHistogram(spike_plot=self)
        self.set_sth_ratio(self._sth_ratio) # Reset
        self._sum = spikesum.SpikeSum(spike_plot=self)
        self.set_sum_ratio(self._sum_ratio) # Reset
        self._axes_pad = 0.025
        self._data_pad = 0.01
        self._draw = True
        
    def plot_spikes(self, spikes, label=None, draw=True, savefig=None, \
            cell_offset=0):
        """
        Plot the spikes. This sets or appends the spike data to be drawn. If
        a single snapshot is needed, simply pass in the desired spikes to plot. 
        If, however, the plot is to be interactive, it is best to send the 
        full range of spikes and then zoom and translate by modifying the 
        axes via the :func:`update_xlim` method.
        
        :param spikes: a 2D vector where cells are in the first dimension and
            their spike times are in the second. This can be a list of lists
            or a 2D numpy.ndarray. If the spike train of a single cell
            is desired for plotting, simply encapsulate the 1D
            vector in square brackets: ``sp.plot_spikes([spikevec])``.
            
        :param label: is the name assigned to this set of spikes. This label
            is used largely internally by the SpikePlot object to 
            differentiate multiple sets of spikes. Default is ``None``. Any
            previously assigned spikes with this label will be replaced.
            
            .. note:: ``None`` is still a label, so repeated calls to 
                :func:`plot_spikes` without explicitly setting the label will
                replace any previously set spikes.
            
        :param draw: is a boolean to specify whether or not to draw. By 
            default, this is ``True``, but it may be desirable to not draw when 
            ascribing multiple sets of spikes and to draw to the screen when 
            setting the last spikes. If other postprocessing is done to the
            figure or axes, it may also make sense to set this value to 
            ``False``.
            
        :param savefig: is a boolean to specify whether to save the figure to
            a file. By default this is ``None`` and is undefined. If, however,
            it is set to ``True`` or ``False``, then this is equivalent to
            calling :func:`set_savefig` with that value.
            
        :param cell_offset: is a vertical offset for this set of spikes in the
            spike raster. By default this is ``0``. However, if the spike 
            raster contains stacked cells of different types (i.e. multiple 
            spike sets each with a different label), then the offset values 
            specifies the vertical offset for these cells.
        """
        # Make sure this is a valid set of spikes
        if isinstance(spikes, list) is False and \
                isinstance(spikes, numpy.ndarray) is False:
            errstr = 'spikes needs to be a list of lists or 2d numpy array'
            raise TypeError(errstr)
            
        # If there are no spikes to draw, just return
        if len(spikes)==0:
            return
        
        # Remove any pre-existing lines
        self._remove_lines(label)
        
        # Get the minimum and maximum values of these spikes
        try:
            spikes_min_x, spikes_max_x = spiketrain.get_spike_bounds(spikes)
        except (TypeError, listutil.ListEmptyError) as ex:
            print(ex)
            #raise(ex)
            return
        
        # If refresh is True, then we completely redraw the figure.
        if self._fig is None:
            # Make a new figure
            self._fig = pyplot.figure(figsize=self._figsize)
            self._raster_axes = None
        if self._raster_axes is None:
            # Make the spike axes. Set to fill the figure and adjust later if
            # other axes are drawn.
            self._raster_axes = self._fig.add_subplot(111)
            self._raster_axes.tick_params(direction='out', length=4, 
                                          width=.75, color='black')
            
        # Update the spike_params dictionary with these spikes.
        try:
            self._spike_params.pop(label)
        except KeyError:
            pass # Ignore if not present
            
        if cell_offset > 0:
            new_height = len(spikes) + cell_offset
            for spike_params in self._spike_params.values():
                msize = self._calculate_markersize(new_height, 
                                                   spike_params.markerscale)
                for line in spike_params.raster_lines:
                    line.set_markersize(msize)

        spike_params = spikeparams.SpikeParams(spikes=spikes, label=label, \
                data_xlim = [spikes_min_x, spikes_max_x], \
                cell_offset=cell_offset, \
                marker=self._marker, markercolor=self._markercolor, \
                markerscale=self._markerscale, \
                markeredgewidth=self._markeredgewidth, \
                linestyle=self._linestyle, linewidth=self._linewidth, \
                sth_redraw=self._sth.redraw, \
                sth_style=self._sth._style, sth_dt=self._sth._dt, \
                sth_kernel=self._sth._kernel, sth_origin=self._sth._origin, \
                sum_linewidth=self._sum._linewidth, \
                sum_redraw=self._sum.redraw, \
                sum_style=self._sum._style)
        self._spike_params[label] = spike_params
        
        # Update the complete data bounds
        self._calculate_data_xlim()

        # Update raster_axes ylim
        height = self._calculate_ylim()
        
        # Get the x data length
        range_x = self._data_xlim[1] - self._data_xlim[0] + 1
        # Pad the visible axis range a bit.
        self._axes_lim[0] = numpy.rint(self._data_xlim[0] - \
                (range_x*self._data_pad))
        self._axes_lim[1] = max( \
        numpy.rint(self._data_xlim[1] + (range_x*self._data_pad)), \
                self._data_xlim[1] + 1)
        self._raster_axes.set_xlim(self._axes_lim)
     
        # markersize to be spike axes height / (num_cells+1)
        m_size = self._calculate_markersize(height)

        raster_lines = []
        # Actually draw the spikes.
        for i in range(len(spikes)):
            row = numpy.ones_like(spikes[i])*(i + 1 + cell_offset)
            line = lines.Line2D(spikes[i], row, \
                    linestyle=self._linestyle, linewidth=self._linewidth, \
                    color=self._markercolor, \
                    marker=self._marker, markeredgecolor=self._markercolor, \
                    markerfacecolor=self._markercolor, markersize=m_size, \
                    markeredgewidth=self._markeredgewidth)
            self._raster_axes.add_line(line)
            raster_lines.append(line)
    
        self._spike_params[label].raster_lines = raster_lines

        # If there is a SpikeTimeHistogram object, create or update it.
        if self._sth_ratio > 0.:
            if self._sth._axes is None:
                self.set_sth_ratio(self._sth_ratio)
            self._sth.plot(spike_params)

        # If there is a SpikeSum object, create or update it.
        if self._sum_ratio > 0.:
            if self._sum._axes is None:
                self.set_sum_ratio(self._sum_ratio)
            self._sum.plot(spike_params)
                    
        # If we are supposed to draw, go ahead.
        self._do_draw(draw)
            
        # If we are supposed to save the figure, do so.
        self._do_savefig(savefig)
    
    def _remove_lines(self, label):
        """Remove the drawn spikes with this label from the raster axes."""
        if self._raster_axes is not None:
            try:
                spike_params = self._spike_params.pop(label)
                for line in spike_params.raster_lines:
                    self._raster_axes.lines.remove(line)
            except KeyError:
                pass # Ignore if the lines do not already exist
                
    def _calculate_data_xlim(self):
        """Calculate the earliest and latest spike times out of all spikes
        in the spikes dictionary."""
        min_x = sys.maxsize
        max_x = 0.
        
        for spike_param in self._spike_params.values():
            xlim = spike_param.data_xlim
            #print xlim
            if xlim[0] < min_x:
                min_x = xlim[0]
            if xlim[1] > max_x:
                max_x = xlim[1]
        self._data_xlim = [min_x, max_x]

    def _calculate_ylim(self):
        """Calculate the total number of rows to draw. """
        total_cells = 0
        for spike_param in self._spike_params.values():
            num_cells = spike_param.cell_offset + len(spike_param.spikes)
            if num_cells > total_cells:
                total_cells = num_cells
        ylim = float(total_cells + 1)
        self._raster_axes.set_ylim(0.5, ylim-.5)
        return total_cells
                
    def _do_draw(self, draw=None):
        """Actually draw the figure. """
        if draw is not None and (draw is True or draw is False):
            self._draw = draw
        if self._draw:
            pyplot.draw()
            pyplot.show()

    def _do_savefig(self, savefig=None):
        """Actually save the figure to a file. """
        if savefig is not None and vartest.isbool(savefig, 'savefig'):
            self._savefig = savefig
        if self._savefig:
            if self._fig_name.endswith('.pdf'):
                pdf = PdfPages(self._fig_name)
                pdf.savefig(self._fig)
                pdf.close()
            else:
                pyplot.savefig(self._fig_name)
                
    def update_xlim(self, xlim):
        """Update the xlim of the raster axes and propagate this info to
        the other axes if they are visible, redrawing and re-writing a
        figure if necessary. 
        
        :param xlim: is a tuple of length 2 specifying the new axes bounds.
        """
        if type(xlim) is tuple and len(xlim)==2:
            self._raster_axes.set_xlim(xlim)
            if self._sth is not None:
                self._sth.update_xlim()
            if self._sum is not None:
                self._sum.update_xlim()
            self._do_draw()    # Draw if necessary
            self._do_savefig() # Save the figure, if necessary
                        
    def _calculate_markersize(self, num_cells, markerscale=None):
        """Given the markerscale and the size of the raster axes, determine
        the size of the marker to draw.
        """
        pos = self._raster_axes.get_position()
        if markerscale is None:
            markerscale = self._markerscale
        return float(self._raster_axes.bbox.size[1])/ \
            float(num_cells) * (pos.height + pos.ymin) * markerscale

    def get_savefig(self):
        """Returns the flag of whether a file is saved upon plotting."""
        return self._savefig
        
    def get_marker(self):
        """Returns the marker type for raster tick marks."""
        return self._marker
    
    def get_markercolor(self):
        """Returns the marker and line color."""
        return self._markercolor
        
    def get_markerscale(self):
        """Returns the relative size of the marker in the raster plot."""
        return self._markerscale
        
    def get_markeredgewidth(self):
        """Returns the size of the tick mark edge. In the case of a vertical
        bar, this is the linewidth of that bar."""
        return self._markeredgewidth
        
    def get_linestyle(self):
        """Returns the linestyle of the raster plot."""
        return self._linestyle
    
    def get_linewidth(self):
        """Returns the line width of the raster plot if the linestyle is not
        'none'."""
        return self._linewidth
        
    def get_fig_name(self):
        """Returns the name of the file that will be written when plotted 
        with the *savefig* variable set to True."""
        return self._fig_name
        
    def get_fig(self):
        """Returns a handle to the Figure."""
        return self._fig
        
    def get_raster_axes(self):
        """Returns a handle to the raster axes."""
        return self._raster_axes
        
    def get_figsize(self):
        """Returns the size of the figure in inches."""
        return self._figsize
        
    def get_sth(self):
        """Returns the instance of the SpikeTimeHistogram object."""
        return self._sth
        
    def get_sum(self):
        """Returns the instance of the SpikeSum object."""
        return self._sum
        
    def get_axes_pad(self):
        """Returns the amount of padding between the raster axes and the
        spike time histogram and/or spike sum axes if they are drawn."""
        return self._axes_pad
        
    def set_axes_pad(self, pad):
        """
        Set the amount of padding between the spike axes and
        the spike time histogram and cumulative spike plot if
        they are shown. 
        
        :param pad: amount to pad. This value is in figure 
            coordinates and needs to be between 0 and 1, and 
            probably very close to 0. A value of 0 means that
            the axes touch each other.
        
        :type pad: float; default 0.01
        """
        vartest.inrange(pad, 0., 1., 'pad')
        self._data_pad = pad

    def get_data_pad(self):
        """Returns the amount of padding in terms of the data width
        to draw a full-screen image. The data width is the last
        spike in the data minus the first spike in the data. 
        The pad ensures that the first and last spike will appear
        in the figure and not be blocked by the vertical axis."""
        return self._data_pad
        
    def set_data_pad(self, pad):
        """
        Set the amount of padding in terms of the data width
        to draw a full-screen image. The data width is the last
        spike in the data minus the first spike in the data. 
        The pad ensures that the first and last spike will appear
        in the figure and not be blocked by the vertical axis. 
        
        :param pad: amount to pad. The default value of 0.025 means
            that the maximum x-axis dimensions will be 
            (first_spike - 0.025*data_width, last_spike + 0.025*data_width).
        
        :type pad: float; default 0.025
        """
        vartest.greater_than_or_equal(pad, 0., 'pad')
        self._axes_pad = pad
                
    def set_marker(self, marker):
        """
        Sets the tick mark symbol. Valid markers are defined by
        `Matplotlib.lines.Line2D markers \
        <http://matplotlib.sourceforge.net/
        api/artist_api.html#matplotlib.lines.Line2D.set_marker>`_.
        """
        try:
            self._line2d.set_marker(marker)
            self._marker = marker
        except Exception:
            print("Invalid marker parameter:", marker, \
            " Use Matplotlib.lines markers.")
            errstr = '%s%s' % ('See http://matplotlib.sourceforge.net/api/', \
                    'artist_api.html#matplotlib.lines.Line2D.set_marker\n')
            print(errstr)
            raise
            
    def set_markercolor(self, color):
        """
        Set the marker color.
        
        :param color: Any `Matplotlib color \
            <http://matplotlib.sourceforge.net/api/colors_api.html>`_.
        
        .. note:: Any subsequent drawing will keep this color,
            so call :func:`get_markercolor` if you want to restore drawing.
        """
        try:
            colorConverter.to_rgba_array(color)
            self._markercolor = color
        except Exception:
            print("Invalid color parameter:", color, " Use Matplotlib colors.")
            print("See http://matplotlib.sourceforge.net/api/colors_api.html\n")
            raise
    
    def set_markerscale(self, markerscale):
        """
        Sets the scale of the tick marks between rows.
        
        :param markerscale: the relative marker size between rows.
            A value of 1 means that the tick mark, if the marker is
            `|`, will make it so that each row just touches. Smaller
            values will make the marker size smaller and larger
            values will cause the markers to bleed across rows.
        """
        vartest.greater_than(markerscale, 0., 'markerscale')
        self._markerscale = markerscale
            
    def set_markeredgewidth(self, markeredgewidth):
        """
        Sets the edge thickness of the tick mark.
        
        :param markeredgewidth: the thickness of the marker.
        
        :type markeredgewidth: float; default 0.75
        """
        vartest.greater_than(markeredgewidth, 0., 'markeredgewidth')
        self._markeredgewidth = float(markeredgewidth)

    def set_linestyle(self, linestyle):
        """
        Sets the horizontal line style. Valid styles are defined by
        `Matplotlib.lines.Line2D linestyle \
        <http://matplotlib.sourceforge.net/api/artist_api.html#matplotlib.lines.Line2D.set_linestyle>`_.
        """
        try:
            self._line2d.set_linestyle(linestyle)
            self._linestyle = linestyle
        except Exception:
            print("Invalid linestyle parameter:", linestyle, \
            " Use Matplotlib.lines linestyle.")
            print("See http://matplotlib.sourceforge.net/api/" \
            "artist_api.html#matplotlib.lines.Line2D.set_linestyle\n")
            raise
            
    def set_linewidth(self, linewidth):
        """
        Sets the width of the horizontal line for each cell, if the 
        linestyle is not 'None'.
        
        :param linewidth: the width of the line.
        
        :type linewidth: float; default 0.75
        """
        vartest.greater_than(linewidth, 0., 'linewidth')
        self._linewidth = linewidth
        
    def set_savefig(self, savefig):
        """
        Set the flag to save a figure after drawing. The name of the file is
        defined by the ``fig_name`` variable.
        
        :param savefig: is a boolean to specify whether to save the plotted 
            figure to a file. By default this is *False*. When set to *True*, 
            the file specified by the ``fig_name`` member variable is written.
    
        .. note:: Since the interactive reference running on the 
            `NEURON SAGE server <https:sage.med.yale.edu:8000>`_
            draws graphics by importing images into the web page. Many 
            examples set ``savefig`` to *True* during instanciation for 
            proper functionality with that guide.

        """
        vartest.isbool(savefig, 'savefig')
        self._savefig = savefig
        
    def set_fig_name(self, fig_name="spikeplot.png"):
        """
        Sets the name of the output figure.
        
        :param fig_name: the output file name.
        
        :type fig_name: str; default "spikeplot.png"
        
        .. note:: The extension should be included with the name and should
            be a valid graphics format supported by the backend drawing
            system. This is not checked for, but most graphics formats
            include png, pdf, ps, eps and svg.
        """
        if type(fig_name) is not str:
            errstr = 'fig_name is %(a)s and of type %(b)s.\n' % \
            {'a':str(fig_name), 'b':type(fig_name)}
            errstr += '  Instead, fig_name should be of type \'str\'\n'
            raise TypeError(errstr)
        self._fig_name = fig_name
        
    def set_fig(self, fig=None):
        """
        Set the figure handle.
        
        :param fig: the figure handle. If this is not None, then this 
            also instantiates the `ax` axes by either creating one if  
            it does not already exist in the figure, or by assigning 
            the first axes object. If None, then the figure handle and
            the axes handle are both set to None.
        """
        if fig is not None:
            if type(fig) is not pyplot.Figure:
                errstr = 'fig is of type %(a)s and should be of type \' \
                matplotlib.figure.Figure\'\n' % \
                {'a':type(fig) }
                raise TypeError(errstr)
            self._fig = fig
            if len(self._fig.get_axes()) == 0 and self._raster_axes is None:
                self._raster_axes = self._fig.add_subplot(111)
            else:
                if self._raster_axes is not None:
                    self._fig.add_axes(self._raster_axes)
                else:
                    self._raster_axes = self._fig.get_axes()[0]
        else:
            self._fig = None
            self._raster_axes = None
            
    def set_raster_axes(self, axes):
        """Sets (or resets) the raster axes to the axes handle passed in.
        This may be necessary when resizing the image or axes so that the
        markers can be resized to scale with the image.
        """
        self._raster_axes = axes
        if axes is not None:
            m_size = self._calculate_markersize(len(axes.get_lines()))
            for line in self._raster_axes.get_lines():
                line.set_markersize(m_size)
        
            # Make sure that the SpikeTimeHistogram has the same width as
            # we do.
            pos = axes.get_position()
            if getattr(self, '_sth', None) and self._sth is not None and \
                    self._sth._axes is not None:
                sth_pos = self._sth._axes.get_position()
                self._sth._axes.set_position([pos.xmin, sth_pos.ymin, \
                        pos.width, sth_pos.height])
                   
            # Make sure that the SpikeSum histogram has the same height as
            # we do.
            if getattr(self, '_sum', None) and self._sum is not None and \
                    self._sum._axes is not None:
                sum_pos = self._sum._axes.get_position()
                self._sum._axes.set_position([sum_pos.xmin, pos.ymin, \
                        sum_pos.width, pos.height])
            
    def set_figsize(self, figsize=None):
        """
        Set the figure size in inches.
        
        :param figsize: a tuple specified as (width, height). If
            None, then it will use the Matplotlib default for a figure
            size when creating a new figure.
        
        :type figsize: tuple; default None
        """
        try:
            pyplot.Figure(figsize) # See if this generates an error.
            self._figsize = figsize
            # Update the markersizes
        except Exception:
            errstr = 'figsize == %(a)s and is of type %(b)s.\n' % \
            {'a':str(figsize), 'b':type(figsize)}
            errstr += '  Instead, figsize should be a \'tuple\' of length 2.\n'
            print(errstr)
            self._figsize = None
            raise
            
    def set_sth_ratio(self, sth_ratio):
        """
        Set the ratio of the spike time histogram axes to the figure.
        
        :param sth_ratio: is the vertical amount between [0,1]
            in figure space to draw a spike time histogram. A value
            of zero means that the spike time histogram axes is not
            drawn. A value of 0.5 means that the figure will split
            vertically in half the spike raster and the histogram.
            A value of 1.0 means that the spike raster is not drawn
            but the histogram takes up the whole space.
            
        :type sth_ratio: float; default is zero.
        """
        vartest.inrange(sth_ratio, 0., 1., 'sth_ratio')
        self._sth_ratio = sth_ratio
        
        if getattr(self, '_fig', None) is None or self._fig is None:
            return # Wait until later to initialize
            
        if sth_ratio > 0.:
            if self._raster_axes is None:
                return

            # Get the raster_axes' current bounds
            pos = self._raster_axes.get_position()

            # Move the bottom up, and shrink the height
            fsize = self._fig.get_size_inches()
            fig_ratio = fsize[0] / fsize[1]
            pad = self._axes_pad
            if fig_ratio > 1.:
                pad *= fig_ratio
            self._raster_axes.set_position([pos.xmin, \
            pos.ymin + (self._sth_ratio * pos.height) + pad, \
            pos.width, \
            pos.height - (self._sth_ratio * pos.height) - pad])
    
            self._sth.set_axes(self._fig.add_axes( \
            [pos.xmin, pos.ymin, pos.width, self._sth_ratio * pos.height]))
            self.set_raster_axes(self._raster_axes)
        
        else:
            # Delete the SpikeTimeHistogram object and axes.
            if self._sth is None:
                return # Nothing to do. Already does not exist
                
            # Resize the raster_axes
            if self._sth._axes is not None:
                pos_sth = self._sth._axes.get_position()
                pos = self._raster_axes.get_position()
                self._raster_axes.set_position([pos.xmin, pos_sth.ymin, \
                pos.width, pos.height + pos_sth.height])
            
            # Turn on the x-axis tick labels
            if self._raster_axes is not None:
                labels = self._raster_axes.get_xticklabels()
                for item in labels:
                    item.set_visible(True)

            # Delete
            if self._fig is not None and self._sth is not None and \
            self._sth._axes is not None:
                self._fig.delaxes(self._sth._axes)
                self._sth_axes = None # Necessary?
 
    def set_sum_ratio(self, sum_ratio):
        """
        Set the ratio of the spike sum histogram axes to the figure.
        
        :param sum_ratio: is the horizontal amount between [0,1]
            in figure space to draw a spike time histogram. A value
            of zero means that the spike time histogram axes is not
            drawn. A value of 0.5 means that the figure will split
            vertically in half the spike raster and the histogram.
            A value of 1.0 means that the spike raster is not drawn
            but the histogram takes up the whole space.
            
        :type sum_ratio: float; default is zero.
        """
        vartest.inrange(sum_ratio, 0., 1., 'sum_ratio')
        self._sum_ratio = sum_ratio
        
        if getattr(self, '_fig', None) is None or self._fig is None:
            return # Wait until later to initialize
            
        if sum_ratio > 0.:
            if self._raster_axes is None:
                return

            # Get the raster_axes' current bounds
            pos = self._raster_axes.get_position()
            fsize = self._fig.get_size_inches()
            fig_ratio = fsize[0] / fsize[1]
            pad = self._axes_pad
            if fig_ratio < 1.:
                pad *= fig_ratio
                
            # Move the right in, and shrink the width
            new_width = pos.width - (self._sum_ratio * pos.width) - pad
            sum_width = self._sum_ratio * pos.width
            self._raster_axes.set_position([pos.xmin, pos.ymin, \
            new_width, pos.height])
    
            self._sum.set_axes(self._fig.add_axes( \
            [pos.xmin+new_width+pad, pos.ymin, \
            sum_width, pos.height]))
            self.set_raster_axes(self._raster_axes)

        else:
            # Delete the SpikeSum object and axes.
            if self._sum is None:
                return # Nothing to do. Already does not exist
                
            # Resize the raster_axes
            if self._sum._axes is not None:
                pos_sum = self._sum._axes.get_position()
                pos = self._raster_axes.get_position()
                self._raster_axes.set_position([pos.xmin, pos_sum.ymin, \
                pos.width, pos.height + pos_sum.height])
            
            # Turn on the x-axis tick labels
            if self._raster_axes is not None:
                labels = self._raster_axes.get_xticklabels()
                for item in labels:
                    item.set_visible(True)

            # Delete
            if self._fig is not None and self._sum is not None and \
            self._sum._axes is not None:
                self._fig.delaxes(self._sum._axes)
                self._sum_axes = None # Necessary?
