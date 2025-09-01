# -*- coding: utf-8 -*-
"""
Decorator and helper routines.

AUTHORS:

- THOMAS MCTAVISH (2010-12-01): initial version, 0.1
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

def set_overrides_and_defaults(obj, defaults, kwargs):
    """For a class that has ``'set_<var>()'`` accessor routines for its variables, this 
    function should be  called in  __init__ to set any kwarg overrides and then set 
    remaining default parameters by calling the appropriate ``'set_<var>()'`` function. 
    To use, simply place default values in a dict and pass it wih the kwargs.
    
    EXAMPLE:
        
    The following example shows general use.
    
    ::
        
        from neuronpy.util.decorators import set_overrides_and_defaults

        class MyClass(object):
            def __init__(self, **kwargs):
                defaults = {'var1':None, \
                        'var2':'spikeplot.png'}
                set_overrides_and_defaults(self, defaults, kwargs)
                    
            def set_var1(self, val):
                # Do some tests to make sure val is appropriate, and then set
                self._var1 = val
                
            def set_var2(self, val):
                # Do some tests to make sure val is appropriate, and then set
                self._var2 = val
                
        mc = MyClass()
        print "mc._var1 =", mc._var1, "mc._var2=", mc._var2
        
    This will give the following output::
        
        mc._var1 = None mc._var2= spikeplot.png
        
    We can override individual variables::
        
        mc = MyClass(var2='case2')
        print "mc._var1 =", mc._var1, "mc._var2=", mc._var2
        
    and that will give::
        
        mc._var1 = None mc._var2= case2
        
    Finally, if the user tries to pass an invalid argument, an error is thrown.::
        
        mc = MyClass(not_in_dict='My Invalid Keyword Arg')
        print "mc._var1 =", mc._var1, "mc._var2=", mc._var2
        
    This is the output.::
        
        AttributeError: Unknown property not_in_dict
        
    .. Note::
        
        This function removes the overriden items from the ``defaults`` dict. 
        If the class needs to retain this dict, it should send a copy to 
        this function.
    """
    # Set the value of the parameters passed in.
    for key, val in kwargs.items():
        func = getattr(obj, 'set_' + key, None)
        if func is None or not callable(func):
            raise AttributeError('Unknown property %s'%key)
        func(val)
        defaults.pop(key)
    
    # Set any parameters not passed in to their default values.
    for key, val in defaults.items():
        func = getattr(obj, 'set_' + key, None)
        if func is None or not callable(func):
            raise AttributeError('Unknown property %s'%key)
        func(val)
