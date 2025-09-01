# -*- coding: utf-8 -*-
"""
Access and test hoc variables.

AUTHORS:

- THOMAS MCTAVISH (2010-04-05): initial version

EXAMPLES:

Retrieve the variables associated with an object::

    >>>vars = get_variables(hoc, obj)

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

def get_variables(hoc, obj, vartype=int(0)):
    """Get the variables associated with an object.
    
    :param hoc: is a hoc object such as would be obtained
            when importing neuron.
            
    :param obj: is the object from which to retrieve
            its variable names. This can be a Section, Segment,
            Mechanism, or Point Process
            
    :param vartype: If vartype = 1, 2, or 3, the storage is for 
            PARAMETER, ASSIGNED, or STATE variables respectively. 
            If vartype = 0, the storage is for all three types.

            If vartype = -1, the count and names (and array size) 
            of the GLOBAL variables are accessible, but any other 
            method will generate an error message. 
    http://www.neuron.yale.edu/neuron/static/docs/help/neuron/neuron/classes/mechstan.html
    
    :returns: a dict object where the keys are the names of the
        variables and the value is filled with the value.
    """
    # See if this is a Section or a Mechanism
    if str(type(obj)) == "<type 'nrn.Segment'>":
        return get_variables_from_segment(hoc, obj, vartype)
    elif str(type(obj)) == "<type 'nrn.Section'>":
        return get_variables_from_section(hoc, obj, vartype)
    elif str(type(obj)) == "<type 'nrn.Mechanism'>":
        return get_variables_from_mechanism(hoc, obj, vartype)
    else:
        return get_variables_from_object(hoc, obj, vartype)

    return None

def get_variables_from_segment(hoc, obj, vartype=int(0)):
    """Get the variables associated with a segment.
    
    :param hoc: is a hoc object such as would be obtained
            when importing neuron.
            
    :param obj: is the Segment object from 
            which to retrieve its variable names.
            
    :param vartype: If vartype = 1, 2, or 3, the storage is for 
            PARAMETER, ASSIGNED, or STATE variables respectively. 
            If vartype = 0, the storage is for all three types.

            If vartype = -1, the count and names (and array size) 
            of the GLOBAL variables are accessible, but any other 
            method will generate an error message. 
    http://www.neuron.yale.edu/neuron/static/docs/help/neuron/neuron/classes/mechstan.html
    
    :returns: a dict object where the keys are the names of the
        variables and the value is filled with the value.
    """
    varnames = {}
    varnames['diam'] = obj.diam
    
    for mech in obj:
        varnames[mech.name()] = mech

    return varnames

def get_variables_from_section(hoc, obj, vartype=int(0)):
    """Get the variables associated with a segment.
    
    :param hoc: is a hoc object such as would be obtained
            when importing neuron.
            
    :param obj: is the Section object from 
            which to retrieve its variable names.
            
    :param vartype: If vartype = 1, 2, or 3, the storage is for 
            PARAMETER, ASSIGNED, or STATE variables respectively. 
            If vartype = 0, the storage is for all three types.

            If vartype = -1, the count and names (and array size) 
            of the GLOBAL variables are accessible, but any other 
            method will generate an error message. 
    http://www.neuron.yale.edu/neuron/static/docs/help/neuron/neuron/classes/mechstan.html
    
    :returns: a dict object where the keys are the names of the
        variables and the value is filled with the value.
    """
    varnames = {}
    varnames['nseg'] = obj.nseg
    varnames['L'] = obj.L
    varnames['Ra'] = obj.Ra
    
    varnames.update(get_variables_from_segment(hoc, obj(0.5), vartype))
    
    return varnames

def get_variables_from_object(hoc, obj, vartype=int(0)):
    """Get the variables associated with a hoc object.
    
    :param hoc: is a hoc object such as would be obtained
            when importing neuron.
            
    :param obj: is the HocObject object from 
            which to retrieve its variable names.
            
    :param vartype: If vartype = 1, 2, or 3, the storage is for 
            PARAMETER, ASSIGNED, or STATE variables respectively. 
            If vartype = 0, the storage is for all three types.

            If vartype = -1, the count and names (and array size) 
            of the GLOBAL variables are accessible, but any other 
            method will generate an error message. 
    http://www.neuron.yale.edu/neuron/static/docs/help/neuron/neuron/classes/mechstan.html
    
    :returns: a dict object where the keys are the names of the
        variables and the value is filled with the value.
    """    
    obj_name = obj.hname()
    mech_name = obj_name.split('[')[0]

    varnames = {}
    
    hoc('objref ms')
    hoc('strdef ms_var_name')
    
    exec_str = 'ms = new MechanismStandard("%s", %d)' % (mech_name, vartype)
    hoc(exec_str)
    for j in range(int(hoc.ms.count())):
        exec_str = 'ms.name(ms_var_name, %d)' % (j)
        hoc(exec_str)
        varnames[hoc.ms_var_name] = \
                eval("obj.%s" % (hoc.ms_var_name))

    return varnames
    
def get_variables_from_mechanism(hoc, obj, vartype=int(0)):
    """Get the variables associated with a hoc object.
    
    :param hoc: is a hoc object such as would be obtained
            when importing neuron.
            
    :param obj: is the Mechanism or Point Process object from 
            which to retrieve its variable names.
            
    :param vartype: If vartype = 1, 2, or 3, the storage is for 
            PARAMETER, ASSIGNED, or STATE variables respectively. 
            If vartype = 0, the storage is for all three types.

            If vartype = -1, the count and names (and array size) 
            of the GLOBAL variables are accessible, but any other 
            method will generate an error message. 
    http://www.neuron.yale.edu/neuron/static/docs/help/neuron/neuron/classes/mechstan.html
    
    :returns: a dict object where the keys are the names of the
        variables and the value is filled with the value.
    """    
    obj_name = obj.name()
    mech_name = obj_name

    varnames = {}
    
    hoc('objref ms')
    hoc('strdef ms_var_name')
    
    exec_str = 'ms = new MechanismStandard("%s", %d)' % (mech_name, vartype)
    hoc(exec_str)
    for j in range(int(hoc.ms.count())):
        exec_str = 'ms.name(ms_var_name, %d)' % (j)
        hoc(exec_str)
        var_name = hoc.ms_var_name
        # remove _obj_name from the variable name
        var_name = var_name.rpartition('_%s'%obj_name)[0]
        varnames[var_name] = \
                eval("obj.%s" % (var_name))

    return varnames