# -*- coding: utf-8 -*-
"""
Variable testing procedures.

Utility procedures to evaluate variables to determine if they pass simple
tests. These tests are typically to ensure that a variable is of a particular
type and a particular format. These methods raise informative exceptions 
if they fail.

AUTHORS:

- THOMAS MCTAVISH (2010-02-05): initial version

EXAMPLES:

Test that a variable is greater than 0.::

    x = 2
    greater_than(x, 0)
    
This does not raise an exception. This shows that `x` may be an `int` or `x`
may be a `float`. Try a negative value.::

    >>> x = -2
    >>> greater_than(x, 0)
    The variable tested needs to be > 0.
    The variable tested is of <type 'int'> and
    has a wrong value of -2.

In this case, the error reported is fairly informative, but note that it can
be made more informative by telling us the name of the parameter. The next
example tries a different format altogether and passes in the name of the
variable so that the user understands which variable is in error.::

    >>> some_var='abc'
    >>> greater_than(some_var, 0, 'some_var')
    some_var needs to be > 0.
    some_var is of <type 'str'> and has a wrong value of abc.

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

def greater_than(test_var, val=0., name=None):
    """
    Tests to see that `test_var` is a float and is
    greater than a number.
    
    :param test_var: the variable to test
    
    :type test_var: float or int
    
    :param name: the name of the parameter to help inform the user
        which variable is problematic, if it does not pass.
               
    :type name: string; default None
    
    :throws: exception if the test_var does not meet the
        requirements and describes the problem.
    """
    try:
        float(test_var)
        if test_var <= val:
            raise ValueError
    except (TypeError, ValueError, SyntaxError):
        if name is None:
            name = 'The variable tested'
        print('%(a)s needs to be > %f.'%{'a':name, 'b':val})
        print('%(a)s is of %(b)s and has a wrong value of %(c)s.\n' % \
        {'a':name, 'b':type(test_var), 'c':test_var })
        raise
        
def greater_than_or_equal(test_var, val=0., name=None):
    """
    Tests to see that `test_var` is a float and is
    greater than zero.
    
    :param test_var: the variable to test
    
    :type test_var: float or int
    
    :param name: the name of the parameter to help inform the user
        which variable is problematic, if it does not pass.
               
    :type name: string; default None
    
    :throws: exception if the test_var does not meet the
        requirements and describes the problem.
    """
    try:
        float(test_var)
        if test_var < 0:
            raise ValueError
    except (TypeError, ValueError, SyntaxError):
        if name is None:
            name = 'The variable tested'
        print('%(a)s needs to be >= 0.'%{'a':name})
        print('%(a)s is of %(b)s and has a wrong value of %(c)s.\n' % \
        {'a':name, 'b':type(test_var), 'c':test_var })
        raise
        
def inrange(test_var, min=0., max=1., name=None):
    """
    Tests to see that `test_var` is a float and is
    greater than or equal to min and less than or
    equal to max.
    
    :param test_var: the variable to test
    
    :type test_var: float or int
    
    :param min: Minimum value
    :type min: float or int

    :param max: Maximum value
    :type max: float or int

    :param name: the name of the parameter to help inform the user
        which variable is problematic, if it does not pass.
               
    :type name: string; default None
    
    :throws: exception if the test_var does not meet the
        requirements and describes the problem.
    """
    try:
        float(test_var)
        if test_var < min or test_var > max:
            raise ValueError
    except (TypeError, ValueError, SyntaxError):
        if name is None:
            name = ''
        else:
            name = " '" + name + "'"
        print('The variable%(a)s needs to be %(b)f <= x <= %(c)f' % \
                {'a':name, 'b':min, 'c':max})
        print('%(a)s is of %(b)s and has a wrong value of %(c)s.\n' % \
                {'a':name, 'b':type(test_var), 'c':test_var })
        raise
        
        
def isbool(test_var, name=None):
    """
    Tests to see that `test_var` is an explicit boolean (True or False) value.
    
    :param test_var: the variable to test
    
    :type test_var: bool
        
    :param name: the name of the parameter to help inform the user
        which variable is problematic, if it does not pass.
               
    :type name: string; default None
    
    :throws: exception if the test_var does not meet the
        requirements and describes the problem.
    """
    try:
        if not(test_var is True or test_var is False):
            raise ValueError
    except (TypeError, ValueError, SyntaxError):
        if name is None:
            name = ''
        else:
            name = " '" + name + "'"
        print('The variable%s needs to be a True or False value.\n' % (name))
        raise

def isint(test_var, name=None):
    """
    Tests to see that `test_var` is an integer.
    
    :param test_var: the variable to test
    
    :type test_var: int
        
    :param name: the name of the parameter to help inform the user
        which variable is problematic, if it does not pass.
               
    :type name: string; default None
    
    :throws: exception if the test_var does not meet the
        requirements and describes the problem.
    """
    try:
        assert(isinstance(test_var, int))
    except (AssertionError):
        if name is None:
            name = ''
        else:
            name = " '" + name + "'"
        print('The variable%s needs to be an integer.\n' % (name))
        raise

def inlist(test_var, valid_list, name=None):
    if test_var not in valid_list:
        if name is None:
            name = ''
        else:
            name = " '" + name + "'"
        errstr = "The variable%s is equal to %s and needs \
                to be one of %s\n" % (name, test_var, valid_list)
        raise ValueError(errstr)

def is1Dvector(test_var, name=None):
    """
    Tests to see that `test_var` is a 1-dimensional list or numpy array.
    
    :param test_var: the variable to test
            
    :param name: the name of the parameter to help inform the user
        which variable is problematic, if it does not pass.
               
    :type name: string; default None
    
    :throws: exception if the test_var does not meet the
        requirements and describes the problem.
    """
    try:
        if (isinstance(test_var, list) is False and \
                isinstance(test_var, numpy.ndarray) is False) or \
                len(numpy.shape(test_var)) != 1:
            raise TypeError
    except TypeError:
        if name is None:
            name = ''
        else:
            name = " '" + name + "'"
        print('The variable%s needs to be a 1-dimensionl python list or numpy array.\n' \
            % (name))
        raise
