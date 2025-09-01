# -*- coding: utf-8 -*-
"""
Methods for platform information.

@author: - Thomas McTavish
"""
import platform

def get_platform_info():
    """
    Retrieve platform information as a dict.
    
    Code borrowed from the file, ``launch.py`` from the Sumatra package.
    """
    network_name = platform.node()
    bits, linkage = platform.architecture()
    return dict(architecture_bits=bits,
                architecture_linkage=linkage,
                machine=platform.machine(),
                network_name=network_name,
                processor=platform.processor(),
                release=platform.release(),
                system_name=platform.system(),
                version=platform.version())
