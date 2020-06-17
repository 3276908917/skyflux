import flux.rot
import flux.parse
import flux.demo
import flux.generate_model
import flux.stokes

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