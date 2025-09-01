# -*- coding: utf-8 -*-
r"""
Generic NEURON cell object. Provides template abstract methods and a generic 
design for NEURON cell models. 

AUTHORS:

- Thomas McTavish (2010-11-04).
- Thomas McTavish (2011-11-06) Update for cleaner docstring and error reporting.
- Thomas McTavish (2012-02-04) Added _resolve_attr to retrieve attributes, 
    sub-attributes, and sub-elements via strings.
"""

import textwrap
import numpy
from neuron import h as nrn

class Cell(object):
    """
    Generic cell template for NEURON cell objects.
        
    .. Note::
    
    Subclasses can override :func:`__init__` completely, 
    or do some of their own initialization first, and then 
    let :class:`Cell` do its initializing, and then the 
    subclass can finish. For example::
        
        class ChildCell(Cell):
            def __init__(self):
                # Do some stuff
                Cell.__init__(self)
                # Do some more stuff

    You might have to get fancy with arguments::
        
        class ChildCell(Cell):
            def __init__(self, someparameter=someval, **kwargs):
                # Do some stuff
                self.someparameter = someparameter
                Cell.__init__(self, **kwargs)
                # Do some more stuff
    """
    def __init__(self, nameprefix='', **kwargs):
        self.nameprefix = nameprefix
        self.x = 0; self.y = 0; self.z = 0
        self.soma = None
        self.synlist = []
        self.create_sections()
        self.build_topology()
        self.build_subsets()
        self.define_geometry()
        self.define_biophysics()
        self.create_synapses()
        self.gid = 0
    
#    def __getattr__(self, name):
#        return self._resolve_attr(self, name)

    def _resolve_attr(self, obj, attrspec):
        """Permit sub-attribute or sub-element retrieval through dot syntax.
        This allows dynamic programming as strings can be assembled
        to retrieve particular cell subobjects and subelements from the cell.
        For example, to connect a NetStim to a cell object that has a variety 
        of synapses, a method could dynamically construct a 
        string to connect that synapse:
            
            def attach_netstims(cells, loc):
                for cell in cells:
                    stim = nrn.NetStim()
                    # Set stim params
                    nc = nrn.NetCon(stim, getattr(cell, loc))
        
        To attach a NetStim to a cell that has ExpSyn objects on the apical and
        dendritic trees, say, then this could be arranged as:
            
            attach_netstims(cells, 'expsyns.apical_syn')
            attach_netstims(cells, 'expsyns.basal_syn')
            
        Note that the sub elements must have a string label, but single list
        indices can also be addressed. For example ``synlist[2]`` could be
        called with the string "synlist.2".
            
        :param obj: Object to query the attribute
        :param attrspec: String of the attribute or element label to retrieve.
        :return: The attribute or object
        """
        attrssplit = attrspec.split(".")
        attr = attrssplit[0]
        try:
            obj = obj[int(attr)] # In case list element
        except ValueError:
            try:
                obj = obj[attr]
            except (TypeError, KeyError, AttributeError):
                obj = getattr(obj, attr)
        except (TypeError, KeyError, AttributeError):
            obj = getattr(obj, attr)
        if len(attrssplit) > 1:
            attrspec = attrspec.partition(".")[2] # right part of the string.
            return self._resolve_attr(obj, attrspec) # Recurse
        return obj

    def create_sections(self):
        """
        Create the sections of the cell. Remember to do this
        in the form::
            
            self.soma = nrn.Section(name='soma', cell=self)
        """
        errstr = "create_sections() is not implemented.\n"
        errstr += textwrap.dedent(self.create_sections.__doc__)
        raise NotImplementedError(errstr)
    
    def build_topology(self):
        """
        Connect the sections of the cell to build a tree. "0" ends are
        toward the soma and "1" ends are distal. For example, to connect 
        the "0" end of a dendrite to the "1" end of the soma::
            
            self.dend.connect(self.soma(1))  
        """
#        errstr = "build_topology() is not implemented.\n"
#        errstr += textwrap.dedent(self.build_topology.__doc__)
#        raise NotImplementedError(errstr)
        pass # May be a 1-compartment neuron. No need to abstract. 
    
    def define_geometry(self):
        """
        Set the 3D geometry of the cell. The length, diameter, and
        number of segments should be set for each section as well as
        the section's (x,y,z) coordinates, if necessary. For example::
            
            self.soma.L = self.soma.diam = 12.6157 # microns
            self.dend.L = 200                      # microns
            self.dend.diam = 1                     # microns
            self.dend.nseg = 5
        """
        errstr = "define_geometry() is not implemented.\n"
        errstr += textwrap.dedent(self.define_geometry.__doc__)
        raise NotImplementedError(errstr)
    
    def define_biophysics(self):
        """
        Assign the membrane properties across the cell. For example::
            
            for sec in self.all: # 'all' exists in parent object.
                sec.Ra = 100    # Axial resistance in Ohm * cm
                sec.cm = 1      # Membrane capacitance in micro Farads / cm^2
            
            # Insert active Hodgkin-Huxley current in the soma
            self.soma.insert('hh')
            self.soma.gnabar_hh = 0.12  # Sodium conductance in S/cm2
            self.soma.gkbar_hh = 0.036  # Potassium conductance in S/cm2
            self.soma.gl_hh = 0.0003    # Leak conductance in S/cm2
            self.soma.el_hh = -54.3     # Reversal potential in mV
        """
        errstr = "define_biophysics() is not implemented.\n"
        errstr += textwrap.dedent(self.define_biophysics.__doc__)
        raise NotImplementedError(errstr)
    
    def create_synapses(self):
        """
        Create synapses (such as ExpSyn) at various
        segments and add them to self.synlist.
            
        For example, in a ball-and-stick cell with a soma and single
        dendrite, to add an exponentially decaying (tau = 2 ms) synapse 
        in the middle of the dendrite::
            
            syn = nrn.ExpSyn(self.dend(0.5))
            syn.tau = 2
            self.synlist.append(syn)
        """
        pass # Ignore if child does not implement.
    
    def build_subsets(self):
        """
        Build section list iterators. This defines the 'all', 
        SectionList, but subclasses may want to define others. 
        If overriden, call Cell.build_subsets(self) to create the 'all'
        SectionList.
            
        For example, in a subclass with two dendrites, 
        ``self.dend[0]`` and ``self.dend[1]`` already defined, 
        we can add them to a "dendlist" iterator::
            
            Cell.build_subsets(self) # Make 'all' iterator
            self.dendlist = nrn.SectionList()
            self.dendlist.append(self.dend[0])
            self.dendlist.append(self.dend[1])
        """
        self.all = nrn.SectionList()
        self.all.wholetree(sec=self.soma)
    
    def connect2target(self, target, thresh=10):
        """
        Make a new :class:`NetCon` with this cell's membrane
        potential at the soma as the source (i.e. the spike detector)
        onto the target passed in (i.e. a synapse on a cell).
        Subclasses may override with other spike detectors.
        """
        nc = nrn.NetCon(self.soma(1)._ref_v, target, sec = self.soma)
        nc.threshold = thresh
        return nc
    
    def is_art(self):
        """Flag to check if we are an integrate-and-fire artificial cell."""
        return 0
    
    def set_position(self, x, y, z):
        """
        Set the base location in 3D and move all other
        parts of the cell relative to that location.
        """
        for sec in self.all:
            for i in range(int(nrn.n3d())):
                nrn.pt3dchange(i, \
                               x-self.x+nrn.x3d(i), \
                               y-self.y+nrn.y3d(i), \
                               z-self.z+nrn.z3d(i), \
                               nrn.diam3d(i))
        self.x = x; self.y = y; self.z = z
    
    def rotateZ(self, theta):
        """
        Rotate the cell about the Z axis.
        """
        rot_m = numpy.array([[numpy.sin(theta), numpy.cos(theta)], \
                             [numpy.cos(theta), -numpy.sin(theta)]])
        for sec in self.all:
            for i in range(int(nrn.n3d())):
                xy = numpy.dot([nrn.x3d(i), nrn.y3d(i)], rot_m)
                nrn.pt3dchange(i, float(xy[0]), float(xy[1]), nrn.z3d(i), \
                               nrn.diam3d(i))