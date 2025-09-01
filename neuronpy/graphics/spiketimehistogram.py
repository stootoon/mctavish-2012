# -*- coding: utf-8 -*-
"""Spike time histogram of spike train data.

This class works with 
:class:`~neuronpy.graphics.spikeparams.SpikeParams` and
:class:`~neuronpy.graphics.spikeplot.SpikePlot` to place a histogram
of spikes under a raster plot (or by itself if the spike raster plot is
not shown).

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
    sp = spikeplot.SpikePlot(sth_ratio=0.3, savefig=True)
    sp.plot_spikes(spikes)

AUTHORS:

- THOMAS MCTAVISH (2010-03-01): initial version, 0.1

- THOMAS MCTAVISH (2010-12-17): Moved non-graphic functionality to spiketrain.
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
from numpy.testing import assert_approx_equal
from matplotlib import pyplot, ticker
from scipy import ndimage

from neuronpy.util import spiketrain
from neuronpy.util import vartest
from neuronpy.math import kernel
from neuronpy.util.decorators import set_overrides_and_defaults

class SpikeTimeHistogram(object):
    def __init__(self, spike_plot, **kwargs):
        self._spike_plot = spike_plot # SpikePlot object
        self._drawn_lines = dict() # axes drawn lines
        self._axes = None
        self.redraw = True
        
        defaults = {'dt' : 1, \
                'kernel' : [1], \
                'origin' : 0, \
                'style' : 'bar', \
                'bins' : None }
                
        set_overrides_and_defaults(self, defaults, kwargs)
        
    def set_bins(self, bins):
        self._bins = bins
        
    def set_axes(self, axes):
        """Set the axes of the SpikeTimeHistogram to the handle passed in."""
        self._axes = axes
        
    def plot(self, spike_params):
        """Calculate and plot the histogram from the SpikeParams object
        passed in."""
        if self._axes is None:
            #print "sth axes is None"
            fig = pyplot.figure()
            self._axes = fig.add_subplot(111)
        
        spike_params.flattened = spiketrain.get_flattened(\
                spike_params.spikes)
            
        try:
            for line in self._drawn_lines[spike_params.label]:
                self._axes.lines.remove(line)
        except KeyError:
            pass # Ignore if not present
        
        # Set the labels
#        self._axes.set_yticks([])
#        self._axes.set_frame_on(False)
    
        self._axes.set_xlim(self._spike_plot._raster_axes.get_xlim())
        
        # Update the data
        self.update_xlim()
        
    def update_xlim(self):
        """Update the xlim of the axes and recalculate the histogram."""
        if self._axes is None:
            return
            
        # Clear the axes and redraw everything
        #self._axes.cla()
    
        # Set the x limits and labels.
        xlabels = self._spike_plot._raster_axes.get_xticklabels()
        for item in xlabels:
            item.set_visible(False)
        xlim = self._spike_plot._raster_axes.get_xlim()
        
        # Go through spike_param dict and draw those histograms.
        for spike_param in self._spike_plot._spike_params.itervalues():
            if spike_param.sth_style is not None:
                if spike_param.sth_redraw:
                    spikes = spiketrain.get_spikes(\
                            [spike_param.flattened], window=(xlim[0], xlim[1]))
                    if self._bins is None:
                        bins = numpy.arange(xlim[0], xlim[1]+self._dt/2, self._dt)
                    else:
                        bins = self._bins.copy()
                        
                    (discrete, bin_edges) = numpy.histogram(spikes, bins)
                    self._filtered = hist_filter(numpy.asfarray(discrete), \
                            spike_param.sth_kernel, spike_param.sth_origin)
                            
                    bins = []

                    for i in xrange(len(self._filtered)):
                        bins.append((i * spike_param.sth_dt) + xlim[0] + \
                                (spike_param.sth_dt / 2))
                    spike_len = len(discrete)
                    rlen = xlim[1] - xlim[0]
                    divisor = (spike_len - \
                            ((spike_len * spike_param.sth_dt) % \
                            float(rlen) / spike_param.sth_dt)) * 1.1

                    if self._style is 'bar':
                        divisor *= 1.15
                    lwidth = 0.75
                    if numpy.float32(divisor) > 0.:
                        lwidth = max(self._axes.bbox.width/divisor, lwidth)

                    if self._style is 'bar' or self._style is 'stepfilled':
                        self._drawn_lines[spike_param.label] = \
                                self._axes.vlines(bins, 0, self._filtered, \
                                        linewidth = lwidth)

                    elif self._style is 'step':
                        x_coords = []
                        y_coords = []
                        for i in xrange(len(self._filtered)):
                            xval = (i * spike_param.sth_dt) + xlim[0]
                            if not(i==0):
                                x_coords.append(xval)
                            x_coords.append(xval)
                            y_coords.append(self._filtered[i])
                            y_coords.append(self._filtered[i])
                        xval = (len(self._filtered) * spike_param.sth_dt) + \
                                xlim[0] + (spike_param.sth_dt / 2)
                        x_coords.append(xval)
                        self._drawn_lines[spike_param.label] = \
                                self._axes.plot(x_coords, y_coords, \
                                        linewidth=spike_param.sth_linewidth, \
                                        color=spike_param.markercolor)

                    elif self._style is 'lineto':
                        x_coords = []                        
                        for i in xrange(len(self._filtered)):
                            xval = (i * spike_param.sth_dt) + xlim[0]
                            x_coords.append(xval)
                        self._drawn_lines[spike_param.label] = \
                                self._axes.plot(x_coords, self._filtered, \
                                        linewidth=spike_param.sth_linewidth, \
                                        color=spike_param.markercolor)
                            
        # Turn off the x-axis tick labels, except the first and last one
        if self._axes is not None:
            self._axes.tick_params(direction='out', length=4, width=.75, color='black')
            self._axes.set_xlim(self._spike_plot._raster_axes.get_xlim())
            self._axes.yaxis.set_major_locator(ticker.LinearLocator(numticks=3))
            self._axes.spines['top'].set_color('none')
            self._axes.get_yaxis().set_ticks_position('left')
            self._axes.get_xaxis().set_ticks_position('bottom')
 
    def set_dt(self, dt):
        """Set the time window of each histogram bar.
        :param dt: Width of the bar (integer or float).
        """
        vartest.greater_than(dt, 0., 'dt')
        self._dt = dt
    
    def set_kernel(self, filter_kernel, test_sum_to_one=True, decimals=2):
        """Set the filter kernel. The spike histogram can be further filtered
        with this kernel. This may be useful in situations where the histogram
        has a small time window (dt) and therefore many bins, which then get
        filtered through this kernel. You will probably want to set the style
        to ``'lineto'`` when filtering this way.
        
        EXAMPLE:
        
        Set the kernel to a gaussian.
        
        ::
            
            # Assumes a SpikePlot object, sp, has already been created.
            from neuronpy.math import kernel
            dt = 0.1
            sth = sp.get_sth()
            sth.set_dt(0.1)
            k = kernel.gauss_1d(sigma=2., dt=dt)
            sth.set_kernel(k)
            sth.set_style('lineto')
            
        This example places
        
        :param kernel: is a 1d vector that should sum to 1.0 to keep the overall
            magnitude of the histogram the same.
        :param test_sum_to_one: is a boolean to test that the kernel indeed
            sums to 1.0. The default is True.
        :param decimals: is the resolution of the kernel summing to 1. Default
            is 2, meaning that the sum has to be 0.99 < sum < 1.01 to be valid.
            This is ignored if test_sum_to_one is false.
        """
        if test_sum_to_one:
            assert_approx_equal(numpy.sum(filter_kernel), 1.0, decimals)
        self._kernel = filter_kernel
        
    def set_origin(self, origin):
        """Set the origin of the filter.
        
        :param origin: can be a string or an integer.:
            - *integer*
              is an integer to specify the offset of the kernel
              applied in the convolution. This needs to be between 0 and
              the length of the kernel - 1. By default, this is 0 and therefore
              centers the kernel. A negative value shifts the kernel to the 
              right and a positive value shifts the kernel to the left.
            - *string*
              Can be 'left', 'center', or 'right'
              
        """
        if isinstance(origin, str):
            valid = ['left', 'center', 'right']
            vartest.inlist(origin, valid, 'origin')
            lenk = len(self._kernel)
            if origin=='left':
                self._origin = -int(numpy.ceil(lenk/2.0))
            elif origin=='right':
                self._origin = int(numpy.ceil(lenk/2.0)) - 1
            else:
                self._origin = 0
            return
            
        vartest.isint(origin, 'origin')
        self._origin = origin
        
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
            - ``'lineto'``
              Draws a line between the (x=bin loc, y=bin height). The difference
              with 'step' is that 'step' draws a horizontal edge across the bin.
              This one draws a line from bin to bin. This is especially useful
              when filtering.
        """
        valid = ['bar', 'step', 'stepfilled', 'lineto']
        vartest.inlist(style, valid, 'style')
        self._style = style
            
    def filter_gaussian(self, sigma=50, limit=0.01):
        """Set the kernel to be a gaussian kernel.
        
        :param sigma: is the standard deviation of the Gaussian distribution.
        
        :param limit: is the threshold for the height of the tail in the
            distribution. Because a Gaussian is infinite, this limits the 
            kernel's length. This value is in a normalized distribution where
            the peak of the bell curve is at 1.0. The default value of 0.01
            is therefore over 4 standard deviations away.
        """
        self.set_kernel(kernel.gauss_1d(sigma/self._dt, limit=limit))
        
#def quantize(spikes, timestep, window):
#    """
#    Aggregate the spike times to one time vector, quantized by dt, 
#    in the bounded window.
#    
#    :param spikes: is the sorted, 1d vector of spike times.
#    
#    :param timestep: is the timestep or width of the bin in time.
#    
#    :param window: is a tuple of length 2 defining the minimum and
#        maximum times of the window.
#    """
#    # We build a new zero-vector that adds 1 to each timestep
#    # where there is a spike. Start by converting the spikes
#    # into index values. We store these values in a vector,
#    # which we will iterate over.
#    
#    # Subtract the offset so that the first spike is the zeroth idx.
#    idx_vec = numpy.add(spikes, -window[0])
#    
#    # Divide those numbers by dt to make index values.
#    idx_vec = numpy.divide(idx_vec, timestep)
#    
#    # Round those numbers to the nearest integer
#    idx_vec = numpy.floor(idx_vec)
#    
#    # Create a vector of zeros that is the length of the SpikePlot
#    # data divided by our timestep.
#    #print "aggregate() xlim=",self._spike_plot._axes_lim
#    num_discrete = numpy.ceil(
#    float(window[1] - window[0]-1) / float(timestep))
#    #print "num_discrete=",num_discrete
#    discrete_vec = numpy.zeros(num_discrete)
#    
#    # Go through and add 1 at the appropriate index
#    for idx in idx_vec:
#        try:
#            discrete_vec[idx] += 1
#        except IndexError:
#            print "idx=", idx
#            
#    return discrete_vec
  
def hist_filter(iline, kernel=[1], origin=0):
    """Convolve the 1d input with the kernel.
    
    :param iline: is the 1d input signal, the discretized histogram
        vector.
        
    :param kernel: is the filter to be applied. This should be a 1d
        vector that sums to 1.0. By default this is [1], which means
        that the output is a copy of the input.
        
    :param origin: is the offset for the filter. If the filter was
        a linear decay, this may need to be set to either the head or
        tail of the filter (i.e. the kernel's length - 1).
        
    ..note::
        
        The time of the histogram bin that is filtered is the left edge and not
        the center of the bin.
    """
    vartest.isint(origin, 'origin')
    lenk = len(kernel)
    min = -int(numpy.ceil(lenk/2.0))
    max = -min - 1
    vartest.inrange(origin, min, max, 'origin')
    return ndimage.convolve1d(iline, weights=kernel, origin=origin, \
            mode='constant')
