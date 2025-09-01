# -*- coding: utf-8 -*-
"""
Initialization of variables passed in the command line or defined in a 
parameter file.

Executes commands defined in a parameters file, as well as processes subsequent 
Python commands delivered on the command line. The usage from the command line 
is then::
    
    python PYTHONFILE [PARAM_FILES &| cmds]

* ``PYTHONFILE`` is the name of a Python script or file.
* ``PARAM_FILES`` is the name(s) of a Python file that execute and may sets 
        global and/or local parameters.
* ``cmds`` are optional Python commands. These can be used to override values 
    set in the parameter file(s) or add new functionality.

EXAMPLE:

Run the simulation defined in ``main.py``, which contains at least the 
following code.

::

    import sys
    from neuron import nrn # a nrn object
    from neuron import h   # a hoc object
    from neuronpy.util import paraminit
    
    # The file ``params.py`` in this same directory contains a dict, 
    # ``sim_var``, which defines our variables.
    from params import sim_var
    
    def run():
        # Run the simulation
        print "Do some amazing science."
        
    if __name__ == "__main__":
        \"\"\"Run when called from the command line.\"\"\"
        local_dict = {}
        paraminit.parse_args(sys.argv[1:], globals(), local_dict)
        if 'sim_var' in local_dict:
            # Replace current values with those in the local_dict
            replacing_dict = local_dict['sim_var']
            for key, val in replacing_dict.items():
                sim_var[key] = val
        run()  # Run the simulation.

From the command line, execute this file with::
    
    python main.py
    
It is important to point out that a file called "params.py" must
also exist at the same level as main.py, unless it is overridden with another
parameter file as in::
    
    python main.py some_params_file.py
    
To pass in other commands, you could deliver something like::
    
    python main.py my_debugging_flag=True 'print "RUNNING IN DEBUGGING MODE"'
    
which will override ``my_debugging_flag`` if it is defined in 
``params.py``, or set it as a new global variable if it is not defined,
and will print the statement before executing the simulation.

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

__version__ = 0.1
import sys, os

def load_param_file(python_file, global_dict, local_dict):
    """
    Load the parameters file(s) into Python's globals dictionary. This is 
    called by :func:`plot_spikes`.
        
    :param python_file: Python file that can be executed with the
            `execfile function 
            <http://docs.python.org/library/functions.html#execfile>`_.
            
    :param global_dict: are the globals from the module that called into this 
        module. Calling ``globals()`` from this module will not carry to the 
        module that is running the simulation.
        
    :param local_dict: is the local dict from the calling module.
        Calling ``locals()`` from this module will not carry to the 
        module that is running the simulation.
    """
    try:
        execfile(python_file, global_dict, local_dict)
    except IOError as (errno, strerror):
        print "I/O error({0}): {1}".format(errno, strerror)
        if errno == 2:
            print "Cannot find the file \'{0}\'".format(python_file)
        raise
    except:
        print "Unexpected error:", sys.exc_info()[0]
        raise

def parse_args(argv, global_dict, local_dict):
    """
    Parse command line arguments and set global parameters. Each argument is
    first tested to see if it is a ``'.py'`` Python file, in which case, it is
    executed. Otherwise, it is executed as a command.
    
    :param argv: is a list of string elements. Since the elements are separated
        by spaces, each element needs to be a contigous Python statement like
        ``variable=value``, with just an equals sign between the variable name 
        and its value and no space. One can use spaces by wrapping an element
        in single quotes as in ``'print "Hello"'``.
        
        If the argument ends with ``'.py'`` it is assumed to be a Python file 
        to execute.ython statement, it is
        assumed to be a filename to a parameter file to load. If no parameters
        file is specified, then the file ``params.py`` in the same 
        directory as the ``PYTHONFILE`` is loaded. If it does not exist, an
        error will be thrown.
    
    :param global_dict: are the globals from the module that called into this 
        module. Calling ``globals()`` from this module will not carry to the 
        module that is running the simulation.
        
    :param local_dict: is the local dict from the calling module.
        Calling ``locals()`` from this module will not carry to the 
        module that is running the simulation.
        
    :return cmd_files, cmds: the list of files and the list of cmds processed.
    """
    cmd_files = []
    cmds = []
    for item in argv:
        if item:
            if os.path.isfile(item) and item.endswith('.py'):
                try:
                    load_param_file(item, global_dict, local_dict)
                except:
                    raise
                cmd_files.append(item)
            else:
                try:
                    exec(item, global_dict, local_dict)
                    cmds.append(item)
                except:
                    warnstr = "WARNING: Do not know how to process statement \
                    \'{0}\'".format(item)
                    print warnstr
                    pass

    return cmd_files, cmds