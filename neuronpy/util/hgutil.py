# -*- coding: utf-8 -*-
"""
Utility methods for interfacing with an hg repository. In particular, this
is useful for programmatically capturing the current (tip) changeset of the
repository, storing that in a file so that if re-running the simulation, you
can revert to the sources used in that simulation.

EXAMPLE::
    from neuronpy.util import hgutil
    
    hg = hgutil.HGUtil('../') # Executed from src directory. Head up to root.
    try:
        hg.check_status()     # Check that all sources have been committed.
    except hgutil.HGStatusNotEmptyError as ex:
        raise ex
    
    tip = hg.get_tip() # Get the last changeset.

    try:
        # Permit arguments to allow a revert to a previous changeset.
        initialize(sys.argv[1:])
        run() # Run the simulation
        finalize() # Write out simulation info
    except:
        # If there is an error, be sure to revert the sources to the last
        # changeset. This does not do anything if the simulatin
        # did not revert to a different changeset.
        hg.revert_to_changeset(str(tip))    

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

import subprocess
from mercurial import hg, ui
from mercurial.node import hex

class HGStatusNotEmptyError(Exception):
    """
    Alert when a hg repository has a non-empty status, meaning that
    there are modifications or changes to the repository from the stored
    changeset.
    """
    def __init__(self, modified=[], added=[], removed=[], deleted=[]):
        self.modified = modified
        self.added = added
        self.removed = removed
        self.deleted = deleted
        
    def __str__(self):
        errstr = 'Changes have been made since the tip changeset.\n'
        errstr += 'You need to either commit these changes or revert\n'
        errstr += 'or rollback to the last changeset!'
        if len(self.modified) > 0:            
            errstr += '\n  The following files have been modified since the \n'
            errstr += '  last changeset:\n'
            for item in self.modified:
                errstr += '    ' + item + '\n'
        if len(self.added) > 0:
            errstr += '\n  The following files have been added since the \n'
            errstr += '  last changeset:\n'
            for item in self.added:
                errstr += '    ' + item + '\n'
        if len(self.removed) > 0:
            errstr += '\n:  The following files have been removed since the \n'
            errstr += '  last changeset:\n'
            for item in self.removed:
                errstr += '    ' + item + '\n'
        if len(self.deleted) > 0:
            errstr += '\n  The following files have been deleted since the \n'
            errstr += '  last changeset:\n'
            for item in self.deleted:
                errstr += '    ' + item + '\n'

        return errstr

class HGUtil(object):
    """
    Utility methods for an hg repository.
    """
    def __init__(self, dir=None):
        self.dir = dir
        self.repo = hg.repository(ui.ui(), dir)
        
    def check_status(self):
        """
        Calls "hg status" and raises a 
        :class:`~neuronpy.util.hgutil.HGStatusNotEmptyError`, if the changeset
        has been modified, or if a file has been added, removed, or deleted.
        """
        modified, added, removed, deleted = self.repo.status()[:4]
        if modified or added or removed or deleted:
            raise HGStatusNotEmptyError(modified, added, removed, deleted)
    
    def get_changeset_dict(self, id=None):
        """
        Return a dictionary of a given changeset given by a changeset id 
        where the keys are the values as returned by the `Mercurial API 
        <http://mercurial.selenic.com/wiki/MercurialApi>_`.
        
        :param id: The changeset ID in hexadecimal format as a string. 
            Default is ``None``, which defaults to ``'tip'``.
        """
        if id is None:
            id = 'tip'
        ctx = self.repo[id]
        
        return dict( \
            rev = ctx.rev(), # the revision number \
            node = hex(ctx.node()), # the revision ID, in hexadecimal \
            user = ctx.user(), # the user who created the changeset \
            date = ctx.date(), # the date of the changeset \
            files = ctx.files(), # the files changed in the changeset \
            description = ctx.description(), #the changeset log message \
            branch = ctx.branch(), # the branch of the changeset \
            tags = ctx.tags(), # a list of the tags applied to the changeset \
            parents = repr(ctx.parents()), # a list of the change context objects  \
                                     # for the changeset's parents \
            children = ctx.children() # a list of the change context \ 
                                      # objects for the changeset's children
            )
            
    def get_tip(self):
        """
        Get the changeset id of the hg tip as a hex number.
        """
        return hex(self.repo.changelog.tip())
        
    def revert_to_changeset(self, id=None):
        """
        Revert the changeset to the id passed in.
        
        :param id: ID of the changeset to revert to. This needs to be a string.
            The default value is ``None``, which reverts to the tip.
            
        :raise subprocess.CalledProcessError: if there is an error.
        """
        try:
            if id is None:
                subprocess.check_call(('hg revert --all %s' % self.dir).split())
            else:
                subprocess.check_call(('hg revert -r %s --all %s' % (id, \
                        self.dir)).split())
        except subprocess.CalledProcessError as ex:
            print ex
            raise ex
