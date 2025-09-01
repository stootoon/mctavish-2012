# -*-(coding: utf-8 -*-
"""
Drawing parameters for a set of spikes used with a 
:class:`~neuronpy.graphics.spikeplot.SpikePlot` object.

AUTHORS:

- THOMAS MCTAVISH (2010-11-01): initial version, 0.1
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

import sys

class SpikeParams(object):
    """
    Each set of spikes that gets drawn gets assigned either default drawing
    properties or those of the SpikePlot instance. They are stored in this
    object.
    """
    def __init__(self, spikes=None, label=None, \
            data_xlim=[sys.maxsize, 0], cell_offset=0, marker=None, \
            markercolor='black', markerscale=1., markeredgewidth=0.75, \
            linestyle='None', linewidth=0.75, sth_redraw=True, \
            sth_style='step', sth_linewidth=0.75, sth_dt=1, \
            sth_kernel=[1], sth_origin=0, sum_redraw=True, sum_style='step', \
            sum_linewidth=0.75):
        self.spikes = spikes
        self.label = label
        self.data_xlim = data_xlim
        self.cell_offset = cell_offset
        self.marker = marker
        self.markercolor = markercolor
        self.markerscale = markerscale
        self.markeredgewidth = markeredgewidth
        self.linestyle = linestyle
        self.linewidth = linewidth
        self.sth_redraw = sth_redraw
        self.sth_style = sth_style
        self.sth_linewidth = sth_linewidth
        self.sth_dt = sth_dt
        self.sth_kernel = sth_kernel
        self.sth_origin = sth_origin
        self.raster_lines = [] # Drawn raster lines in the raster_axes
        self.flattened = []
        self.sum_redraw = sum_redraw
        self.sum_style = sum_style
        self.sum_linewidth = sum_linewidth
