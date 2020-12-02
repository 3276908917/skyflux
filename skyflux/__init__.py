import skyflux.utils

import skyflux.rot
import skyflux.catalog
import skyflux.ant
import skyflux.vis
import skyflux.stokes
import skyflux.demo

"""
I wanted to use this block of code to automatically import any new
scripts that I add, but for some reason it breaks relative imports
in sub-modules.

import pkgutil

__all__ = []

for loader, module_name, is_pkg in pkgutil.walk_packages(__path__):
    __all__.append(module_name)
    _module = loader.find_module(module_name).load_module(module_name)
    globals()[module_name] = _module
"""
