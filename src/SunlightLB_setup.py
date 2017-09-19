#!/usr/bin/python3

# This file is part of SunlightLB - a 3D Lattice Boltzmann code
# Copyright (C) 2005-2014 Unilever UK Central Resources Ltd.

# SunlightLB is free software; you can redistribute it and/or 
# modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation; either version 2 
# of the License, or (at your option) any later version. 

# SunlightLB is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU General Public License for more details. 

# You should have received a copy of the GNU General Public License 
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

"""
SunlightLB_setup.py file for lattice Boltzmann simulations
"""

from distutils.core import setup, Extension

SunlightLB_module = Extension('_SunlightLB', 
                              sources = ['SunlightLB.i'],
                              libraries = ['sunlightlb'],
                              library_dirs = ['.'])

setup(name = 'SunlightLB',
      version = '0.1',
      author = "Patrick Warren",
      description = """Package for lattice Boltzmann simulations""",
      ext_modules = [SunlightLB_module],
      py_modules = ["SunlightLB"])
