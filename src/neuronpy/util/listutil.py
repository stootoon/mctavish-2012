# -*- coding: utf-8 -*-
"""Utility methods for list objects.

AUTHORS:
- Thomas McTavish
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

class ListEmptyError(Exception):
    """Alert when a list is empty"""
    def __init__(self, list_name=None):
        self.list_name = list_name
        
    def __str__(self):
        the_list_name=''
        if self.list_name is not None:
            the_list_name = '\''+self.list_name+'\''
        errstr='List %s is empty.'%the_list_name
        return repr(errstr)
    
def nonempty_copy(list):
    """
    Makes a copy of the list and then removes any empty elements.
    """
    list_copy=list[:]
    remove_empties(list_copy)
    return list_copy

def remove_empties(list):
    """
    Removes any empty elements from the list.
    """
    contains_empties = True
    while contains_empties is True:
        # We may have to re-loop if there are adjacent empty elements.
        contains_empties = False
        for i in iter(list):
            if len(i)==0:
                contains_empties=True
                list.remove(i)
                
def flatten_from_2d(list_of_lists):
    """
    Returns a 1d, flattened version of a 2d array or list of lists.
    This also removes any empty elements.
    """
    vec = []
    for row in iter(list_of_lists):
        vec.extend(row)
        
    return vec