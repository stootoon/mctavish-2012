# -*- coding: utf-8 -*-
"""Spike sum histogram of spike train data.

This class works with 
:class:`~neuronpy.graphics.spikeparams.SpikeParams` and
:class:`~neuronpy.graphics.spikeplot.SpikePlot` to show the cumulative
count of spikes from one cell from raster plots containing rows of 
spike trains.

For detailed examples, see

EXAMPLES:

This example generates some random "spikes" and displays them.

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
    # times. Now, let's plot these spikes with the default parameters.
    sp = spikeplot.SpikePlot(sum_ratio=0.3, savefig=True)
    sp.plot_spikes(spikes)
    
AUTHORS:

- THOMAS MCTAVISH (2010-03-01): initial version, 0.1
"""
# While this software is under the permissive MIT License, 
# (http://www.opensource.org/licenses/mit-license.php)
# We ask that you cite the neuronpy package (or tools used in this package)
# in any publications and contact the author with your referenced publication.
#
# Format:
# McTavish, T.S. NeuronPy library, version 0.1, http://bitbucket.org/tommctavish/neuronpy
#
# Copyright (c) 2010 Thomas S. McTavish
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

import numpy
from matplotlib import pyplot, ticker

from neuronpy.util import spiketrain
from neuronpy.util.decorators import set_overrides_and_defaults
from neuronpy.util import vartest

class SpikeSum(object):
    """Draws a histogram on the right side of a raster plot that is the
    cumulative sum of spikes in the window.
    """
#    def __init__(self, spike_plot):
#        self._spike_plot = spike_plot # SpikePlot object
#        self._axes = None
#        self._style = 'step'
#        self._linewidth = 0.75
#        self._drawn_lines = dict() # axes drawn lines

    def __init__(self, spike_plot, **kwargs):
        self._spike_plot = spike_plot # SpikePlot object
        self._drawn_lines = dict() # axes drawn lines
        self._axes = None
        self.redraw = True
        
        defaults = {'style' : 'step', \
                'linewidth' : 0.75}                
        set_overrides_and_defaults(self, defaults, kwargs)

    def set_axes(self, ax):
        """Set the axes to the axes handle passed in. """
        self._axes = ax
    
    def set_style(self, style):
        """Set the style of the histogram bars.
        :param style: May be one of the following:
            - ``'bar'``
              Makes a solid bar from zero to the bin height (default).
            - ``'stepfilled'``
              Similar to ``'bar'``, but leaves no whitespace between the bars.
            - ``'step'``
              Draws an outline of the bars. This is especially useful if
              drawing multiple histograms.
        """
        valid = ['bar', 'step', 'stepfilled']
        vartest.inlist(style, valid, 'style')
        self._style = style
        
    def set_linewidth(self, linewidth):
        vartest.greater_than_or_equal(linewidth, 0, 'linewidth')
        self._linewidth = linewidth
        
    def plot(self, spike_params):
        """Set the cell spike-count histogram from the SpikeParams object
        passed in.
        
        :param spike_params: is a SpikeParams object. The spikes and line
            drawing parameters contained in it define what is drawn.
        """
        if self._axes is None:
            fig = pyplot.figure()
            self.set_axes(fig.add_subplot(111))


        try:
            self._axes.lines.remove(self._drawn_lines[spike_params.label])
        except KeyError:
            pass # Ignore if not present
        
        # Set the labels
        self._axes.set_yticks([])
        self._axes.set_ylim(self._spike_plot._raster_axes.get_ylim())
        self._axes.spines['right'].set_color('none')
        self._axes.get_yaxis().set_ticks_position('left')

        # Update the data
        self.update_xlim()

    def update_xlim(self):
        """Respond to the fact that the xlimits have changed by redrawing. """
        if self._axes is None:
            return

        # Clear the axes and redraw everything
        self._axes.cla()
        xlim = self._spike_plot._raster_axes.get_xlim()
        #print "xlim",xlim
        
        # Go through spike_param dict and draw those histograms.
        for spike_param in self._spike_plot._spike_params.itervalues():
            if spike_param.sum_style is not None:
                if spike_param.sum_redraw:
                    spikes = spiketrain.get_spikes(spike_param.spikes, \
                            idx=range(len(spike_param.spikes)), \
                            window=(xlim[0], xlim[1]))
                    num_per_train = []
                    for train in spikes.itervalues():
                        num_per_train.append(len(train))
                    if self._style is 'bar' or self._style is 'stepfilled':
    #                    y_coords = range(spike_param.cell_offset,
    #                                     spike_param.cell_offset+len(num_per_train))
    #                    y_coords = numpy.add(y_coords, 1)
    
    
                        ylim = self._spike_plot._raster_axes.get_ylim()
                        lwidth = self._spike_plot._calculate_markersize(ylim[1])
                        if self._style is 'bar':
                            lwidth *= 0.9
                        lwidth = max(lwidth, 0.75)
                        y_coords = []
                        x_coords = []
                        for i in range(len(num_per_train)):
                            x_coords.append(num_per_train[i])
                            yval = ylim[0] + i + spike_param.cell_offset + .5
                            y_coords.append(yval)
                        self._drawn_lines[spike_param.label] = \
                        self._axes.hlines(y_coords, [0], num_per_train, \
                                linewidth = lwidth)
    
                    else: # 'step'
                        ylim = self._spike_plot._raster_axes.get_ylim()
                        y_coords = []
                        x_coords = []
                        for i in range(len(num_per_train)):
                            x_coords.append(num_per_train[i])
                            x_coords.append(num_per_train[i])
                            yval = ylim[0] + i + spike_param.cell_offset
                            if not(i==0):
                                y_coords.append(yval)
                            y_coords.append(yval)
                        y_coords.append(yval+1)
                        self._axes.plot(x_coords, y_coords, \
                                linewidth=spike_param.sum_linewidth, \
                                color=spike_param.markercolor)
                        self._axes.set_xlim(numpy.min(x_coords)-1, numpy.max(x_coords)+1)
                            
            # Turn off the x-axis tick labels, except the first and last one
            if self._axes is not None:
                self._axes.set_ylim(self._spike_plot._raster_axes.get_ylim())
                self._axes.set_yticks([])
                self._axes.tick_params(direction='out', length=4, width=.75, color='black')
                self._axes.xaxis.set_major_locator(ticker.LinearLocator(numticks=3))
            