# -*- coding: utf-8 -*-
r"""
Interface for a simple data store where the root and subtables are nested 
Python :py:class:`dict` objects.

AUTHORS:

- THOMAS MCTAVISH (2010-03-01): initial version, 0.1

EXAMPLES::
    
    from neuronpy.util import dictdb
    import datetime
    
    the_dict = { \
        'sub1' : { \
            'key_a' : { \
                'key_a_a' : 1, \
                'time_stamp' : \
                    datetime.datetime(2010, 7, 23, 18, 43, 36, 640692), \
                'key_a_b' : 'ABCDEFG' }, \
            'key_b' : [1,2,3], \
            'key_c' : { \
                'key_c_a' : 2, \
                'key_c_b' : 'HIJKLMNO' }}, \
        'sub2' : { \
            'key_a' : { \
                'key_a_a' : 2, \
                'time_stamp' : \
                    datetime.datetime(2010, 8, 23, 18, 43, 36, 640692), \
                'key_a_b' : 'XYZPDQ' }, \
            'key_b' : [4,5,6], \
            'key_c' : { \
                'key_c_a' : 1, \
                'key_c_b' : 'ABCDEFG' }}, \
        'different_sub1' : { \
            'different_key_a' : { \
                'different_key_a_a' : 3, \
                'different_key_a_b' : [1,2,3,4,5] }, \
            'different_key_b' : None, \
            'different_key_c' : { \
                'different_key_c_a' : 2.0, \
                'different_key_c_b' : ['a', 'b', 'c'] }} \
    }

Nodes ``'sub1_a'``, ``'sub1_b'``, and ``'different_sub1_a'`` refer to trunk
nodes, which may be thought of as records in the database.
Notice that in this dict, subdictionaries can each be of arbitrary depth, and 
contain any data type that can be put into a dict. This example shows records 
``'sub1_a'`` and  ``'sub1_b'`` share the same structure, but ``different_sub1``
contains different data types.

There are two related functions for retrieving records: :func:`filter_dict` and 
:func:`match`. Both approaches use user-defined :py:keyword:`lambda` functions.
The :func:`filter_dict` function is itself a generator function that operates 
on ``(key, value, parent_keys)`` tuples. This function therefore allows 
any of these values to be ignored and one can search for keys or values 
irrespective of their key associations.

The :func:`match` method permits multiple queries where a given key meets some
condition. The key and the
condition as a :py:keyword:`lambda` function  are provided in a tuple: 
``('key', lambda v: <some_operation_with_v>)`` where *"some_operation_with_v"* 
evaluates to ``True`` or ``False``.
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

def filter_dict(the_dict, predicate=lambda k, v, p: True):
    r"""
    Filter a dict by some function on ``(key, value, parent_key)`` tuples.
    
    :param the_dict: The dict to filter.
    
    :param predicate: A :py:keyword:`lambda` function on ``(k, v, p)``, 
        key, value, key_parent tuples. Where ``k`` and ``v`` are their normal,
        ``key-value``s and where ``p`` is a list of sub-dict keys, from the 
        trunk of the_dict to the current key, ``k``.
        
    EXAMPLES:

    Retrieve a subdict that contains any value with ``'ABCDEFG'``. In this 
    case, ``k`` and ``p`` in the lambda function are ignored.
    
    .. code-block:: python

        filtered = dict(filter_dict(the_dict, lambda k, v, p: v=='ABCDEFG'))
        
    Retrieve a subdict where the ``'time_stamp'`` is greater than 30 days.
    
    .. code-block:: python

        def date_filter(k, v, p):
            return k == 'time_stamp' and \
                    (datetime.datetime.now() - v).days >= 30
                    
        filtered = dict(filter_dict(the_dict, date_filter(k, v, p)))
    
    If a sub-key is not unique, then we can query the parent keys. In this 
    case, ``p`` is a list of keys starting from the first sub-dict to the 
    current key. Therefore, ``k == p[-1]``. This example searches for string
    'ABCDEF', but mandates that ``k == 'key_c_b'`` and the parent key is
    ``p == 'key_c'``.

    .. code-block:: python

        def f(k, v, p):
            return k == 'key_c_b' and len(p) >= 2 and p[-2] == 'key_c'
                    
        filtered = dict(filter_dict(the_dict, f(k, v, p)))

    """
    parent_keys = []
    for k, v in the_dict.iteritems():
        parent_keys[:] = []
        parent_keys.append(k)
        if isinstance(v, dict) and _filter_dict_sub(predicate, v, parent_keys):
            yield k, v
            
def _filter_dict_sub(predicate, the_dict, parent_keys):
    """
    Sub-function to filter_dict for recursing.
    """
    for k, v in the_dict.iteritems():
        parent_keys.append(k)
        if isinstance(v, dict) and _filter_dict_sub(predicate, v, parent_keys):
            parent_keys.pop()
            return True
        if predicate(k, v, parent_keys):
            parent_keys.pop()
            return True
        parent_keys.pop()
    return False
    
def match(the_dict, queries, mode = 'AND', search_type = 0):
    r"""
    Retrieve the root items of ``the_dict`` where 
    nested keys match the conditions to some number of queries.
    
    :param the_dict: The :py:class:`dict` to search.

    :param queries: A tuple or list of 2- (or preferrably 3-) element tuples 
        with the first element being the nested key to search for and the 
        second tuple is the condition on that key's value. If a query cannot 
        be processed, it will be ignored, but a warning will print.
        
        An optional third element can be specified as the list of parent
        keys to help speed the search query and also help guarantee precise 
        results. By default, the first key found in any nested dict that 
        matches the key being searched for in the 
        query (the first element in the tuple) will be the only key evaluated. 
        The search algorithm will not continue looking for the key in other 
        nested dicts. Therefore, duplicate keys in nested dicts should be 
        avoided. If they are used, however, then you can specify the key of 
        the parent.

        .. Note::
            
            It is efficient to specify the parent keys as much as possible
            because the query key may be nested quite deep. Using parent keys
            permits a more direct lookup of the query key.
        
    :param mode: Is ``'AND'`` by default and can be ``'OR'``. To mix,
        ``'AND'`` and ``'OR'`` queries, retrieve the sub-dictionary in one 
        condition and then call this function again passing in that
        sub-dictionary with queries for the other condition.
        
        .. Note::
            
            Two or more ``'AND'`` queries prune a copy of the dictionary. The
            first query returns a subdictionary that the second query operates
            on and so forth. Searches through the subdictionaries will 
            therefore be faster if the first queries prune more of the original
            dict than subsequent queries. Likewise, when mixing ``'AND'`` and 
            ``'OR'`` queries, it is more efficient to process ``'AND'``
            queries first.
            
    :param search_type: By default is 0, but can be 1, or 2. Searches are
        depth-first. This number specifies what is to happen when a key (and
        any parent keys) are found:
            
        - ``0``
            Shallow search. When a match is found exit up to the root dict and 
            do not continue searching the trunk this branch is on.
        
        - ``1``
            Do not continue searching this branch.
            
        - ``2``
            Continue search this branch for sub-keys with the same value.
            

    :return: A sub-dict of ``the_dict`` where nested keys match the conditions
        to queries. The elements in the returned dict are
        copies of ``the_dict``, which remains unaltered.
    
    EXAMPLES:
        
    To search for records where ``'key_a_a' == 2``:
    
    .. code-block:: python
    
        subdict = dictdb.match(the_dict, ('key_a_a',lambda v: v==2))
        for key in subdict.keys():
            print key
    
    This should print ``sub2``.
    
    To make compound queries, pass a list of query-tuples.
        
    .. code-block:: python
    
        subdict = dictdb.match(the_dict, [ \
                ('key_a_a',lambda v: v==2,['key_a']), \
                ('different_key_a_b',lambda v: sum(v)> 3,['different_key_a']) \
                ], \
                mode='OR')
        for key in subdict.keys():
            print key
    
    This should print ``sub2`` and ``different_sub1_a``.
    
    In this case, each query is a list of 3-element tuples. The last
    element in the tuple specifies the parent dict of the particular key being
    searched for. This also makes an OR search of the terms.
    
    For complex or custom data types and expressions, you can define your own 
    function and pass that in as the second parameter in the query. The following
    example tests for a value in a list.
        
    .. code-block:: python
    
        def my_func(x, y):
            return isinstance(x, list) and y in x
        
        subdict = dictdb.match(the_dict, \
                ('different_key_c_b', lambda v: my_func(v, 'b')))
        for key in subdict.keys():
            print key
                
    Here is another example that retrieves records with a timestamp older than
    30 days.
    
    .. code-block:: python
    
        subdict = dictdb.match(the_dict, \
                ('time_stamp', lambda v: (datetime.datetime.now() - x).days >= 30))
        for key in subdict.keys():
            print key
    """
    db_copy = the_dict.copy()
    db_copy2 = {}

    root_matches = []
    calling_keys = []
    sub_dict = {}
    
    if isinstance(queries, tuple):
        queries = [queries]

    for query in queries:
        try:
            if not isinstance(query, tuple):
                raise TypeError
        except (TypeError, ValueError, SyntaxError):
            warning = 'WARNING: Ignoring %(a)s\n' % {'a':query}
            warning += '   It needs to be a 2- or 3-element tuple.'
            print warning
        
        key_str = ''
        if len(query) == 3:
            try:
                if not isinstance(query[2], list) or \
                        not isinstance(query[2][0], str):
                    raise UserWarning
                for key in query[2]:
                    key_str += '[\'' + key + '\']'
            except UserWarning:
                warning = 'WARNING:\n'
                warning += '   \"%(a)s\"\n' % {'a':query[2]}
                warning += '   is an incorrect parent key.\n'
                warning += '   It needs to be a list of string elements\n'
                warning += '   specifying parent keys. Ignoring this list.'
                print warning
                key_str = ''
            except Exception as _ex:
                raise _ex
            
        key_str += '[\'' + query[0] + '\']'
        
        if mode == 'AND': # Erase calling_keys
            calling_keys[:] = []
            root_matches[:] = []

        _search_dict(db_copy, key_str, query[1], root_matches, calling_keys, \
        False, search_type)
        if mode == 'AND':
            db_copy2.clear()
            db_copy2 = {}
            for match_path in root_matches:
                if len(match_path) > 0:
                    the_id = match_path[0]
                    db_copy2[the_id] = db_copy[the_id]                
            db_copy = db_copy2.copy()
            
    for match_path in root_matches:
        if len(match_path) > 0:
            the_id = match_path[0]
            sub_dict[the_id] = db_copy[the_id]
    
    return sub_dict

def _search_dict(the_dict, key_str, condition, root_matches, calling_keys, \
keyfound, search_type):
    """
    Search ``the_dict`` for the given query. Append ``root_matches`` if
    there is a match.
    
    :param the_dict: A dictionary object.
    
    :param key_str: Is a string of keys to attempt to get, which may include
        any number of ancestors. Something in the form
        ``"['grandparent']['parent_key']['key_to_find']"``.
            
    :param root_matches: The cumulative path of keys back to the root, for
        each match is obtained.
    
    :param calling_keys: Record of of the calling keys. If a match is made, it
        gets added to root_matches.
        
    :param keyfound: is a flag to determine if a matching key has been found.
    
    :param search_type: See description in :func:`match`.
    
    :return: True if the key was found, regardless of whether there was a 
        match and False otherwise.
    """        
    try:
        eval_str = 'the_dict' + key_str
        # Test for existence. Will raise a KeyError if it does not exist
        var = eval(eval_str)
        keyfound = True
        
        # Now test the condition
        if condition(var):
            root_matches.append(calling_keys[:])
        if search_type < 2:
            return True # Return True to the fact that we found the key,
                        # not whether the match was made.
    except KeyError:
        pass
    except TypeError:
        raise
    
    # Key not found. Recurse, searching sub-dicts.
    for (key, val) in the_dict.iteritems():
        calling_keys.append(key)
        if isinstance(val, dict):
            try:
                keyfound = _search_dict(val, key_str, condition, root_matches, \
                        calling_keys, keyfound, search_type)
            except TypeError: # Probably tried to explore too deeply
                calling_keys.pop()
                return False
        calling_keys.pop()
        
        if search_type == 0 and keyfound and len(calling_keys)>0:
            return True
        elif search_type == 1:
            keyfound = True
        else:
            keyfound = False
    
    keyfound = False    
    return False
