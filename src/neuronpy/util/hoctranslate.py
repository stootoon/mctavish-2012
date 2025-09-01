# -*- coding: utf-8 -*-
"""
Utility methods for translating Python dicts into hoc template objects or
global variables in NEURON accessible to hoc.

@author: - Thomas McTavish
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

def bool_to_hoc_double(hoc, key, val):
    """Convert a python boolean into a global hoc double value.
    :param hoc: A hoc instance as attained by ``from neuron import h``.
    :param key: Name of the variable
    :param val: True or False
    """
    exec_str = 'hoc(\"' + key + '=' + str(int(val)) + '\")'
    exec(exec_str)

def num_to_hoc_double(hoc, key, val):
    """Convert a python number into a global hoc double value.
    :param hoc: A hoc instance as attained by ``from neuron import h``.
    :param key: Name of the variable
    :param val: Python integer or float
    """
    exec_str = 'hoc(\"' + key + '=' + str(val) + '\")'
    exec(exec_str)
    
def str_to_hoc_strdef(hoc, key, val):
    """Convert a python string into a global hoc strdef.
    :param hoc: A hoc instance as attained by ``from neuron import h``.
    :param key: Name of the variable
    :param val: Python string
    """
    exec_str = 'hoc(\"strdef ' + key + '\")'
    exec(exec_str)
    exec_str = 'hoc.' + key + '=\"' + str(val) + '\"'
    exec(exec_str)
    
def numpy_array_to_hoc_vector(hoc, key, val):
    """Convert a 1D numpy array into a global hoc Vector.
    :param hoc: A hoc instance as attained by ``from neuron import h``.
    :param key: Name of the variable
    :param val: 1D numpy array
    """
    exec_str = 'hoc(\"objref ' + key + '\")'
    exec(exec_str)
    exec_str = 'hoc.' + key + '=hoc.Vector(val)'
    exec(exec_str)

def list_to_hoc_vector(hoc, key, val):
    """Convert a 1D list of ints or floats into a global hoc Vector.
    :param hoc: A hoc instance as attained by ``from neuron import h``.
    :param key: Name of the variable
    :param val: 1D list of ints or floats
    """
    exec_str = 'hoc(\"objref ' + key + '\")'
    exec(exec_str)
    if len(val) > 0:
        floatvec = map(float, val)
        exec_str = 'hoc.' + key + '=hoc.Vector(floatvec)'
        exec(exec_str)

def dict_to_global(hoc, the_dict, show_warning=True):
    """
    Translate numbers, strings, and boolean values of a Python dictionary into 
    global variables in NEURON. The keys of the dict become the variable names 
    in hoc.  If the hoc variables already exist, they will be overwritten. All 
    values that are boolean, ints, or floats are translated to a double in hoc.
    Strings are treated as hoc ``strdef`` objects. Any
    other types of the dict are ignored.
    
    :param hoc: A hoc instance as attained by ``from neuron import h``.
    
    :param the_dict: Dict to translate.
    """
    exlist = []
    for (key, val) in the_dict.iteritems():
        try:
            # Number convert to hoc double
            if type(val) == type(int()) or \
                    type(val) == type(float()):
                num_to_hoc_double(hoc, key, val)
            # Boolean convert to hoc double
            if type(val) == type(bool()):
                bool_to_hoc_double(hoc, key, val)
            # String convert to hoc strdef
            elif type(val) == type(str()):
                str_to_hoc_strdef(hoc, key, val)
            # Numpy array convert to hoc Vector
            elif isinstance(val, numpy.ndarray):
                numpy_array_to_hoc_vector(hoc, key, val)
            # List, convert to hoc Vector
            elif type(val) == type(list()):
                list_to_hoc_vector(hoc, key, val)
        except Exception:
            # Append any that cannot be converted and continue
            exlist.append(key)
            
    if show_warning and len(exlist) > 0:
        warnstr = 'WARNING: In hoctranslate.dict_to_global(), '
        warnstr += 'the following variables were not translated into '
        warnstr += 'hoc globals:\n    ['
        for item in exlist:
            warnstr += item + ', '
        warnstr += ']'
        print warnstr
        
