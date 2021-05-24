#!/usr/local/env python
"""
Generate a shortcut for the application.
"""

import sys
import os

import pymol

# After the installation, create a startmenu shortcut.
# get the python entry file
pymol_file = pymol.__file__

# create a shortcut to the executable file, using the right environment
from pyshortcuts import make_shortcut

if os.name == 'nt':
    # todo: only tested on conda base env
    make_shortcut(sys.executable + ' ' + pymol_file, name='PymolMda',
                  icon='data/pymol/icons/icon2.svg', terminal=False)
else:
    make_shortcut(pymol_file, name='PymolMda',
                  icon='data/pymol/icons/icon2.svg', terminal=False)