import molgif
import os
import random
import ase.build

desk = os.path.join(os.path.expanduser('~'), 'Desktop')
os.chdir(desk)

# load in molecule (ase.Atoms object)
molecule = ase.build.molecule('C60')

# pixels in largest dimension (width or height)
max_px = 900

"""
Use smart_rotate to find the best viewing angle
"""
import ase.io

path = 'C:\\Users\\mcowa\\Box Sync\\Michael_Cowan_PhD_research\\Papers\\02_au30_solubility\\solubility_paper_structures\\au25pet18_ROTATED.xyz'
molecule = ase.io.read(path)

# create rotating gif with rot_gif function
molgif.rot_gif(molecule, os.path.join(desk, 'NEW-AU25'), max_px=max_px, smart_rotate=False,
               save_frames=False, optimize_gif=True, scale=0.8)

import sys
sys.exit()
"""
Add a legend
"""

molgif.rot_gif(molecule, max_px=max_px, add_legend=True)

"""
Specify the color of each atom
"""

# can be a string for one color or a list of custom colors
rainbow = ['red', 'orange', 'yellow',
           'green', 'blue', 'violet'] * 4

# list much match number of atoms
rainbow = rainbow[:len(molecule)]

molgif.rot_gif(molecule, max_px=max_px, colors=rainbow)

"""
Use a dictionary to quickly color by atom type
"""

# default colors will be used for types not specified
molgif.rot_gif(molecule, max_px=max_px, colors=dict(C='hotpink'),
               add_legend=True)

"""
Anchor an atom to be at the center of rotation
"""

# define index of atom to anchor
anchor = 3

colors = ['white'] * len(molecule)
colors[anchor] = '#0892d0'

molgif.rot_gif(molecule, max_px=max_px * 1.54, colors=colors,
               anchor=anchor)

"""
Adjust loop time and FPS
"""

# loop_time = time to complete one rotation (seconds)
molgif.rot_gif(molecule, max_px=max_px, loop_time=2, fps=60)

"""
Turn off bonds and scale atomic sizes
"""

molgif.rot_gif(molecule, max_px=max_px, add_bonds=False,
               scale=0.9)

"""
Change rotation axis
"""

# switch between x, y (Default), or z
molgif.rot_gif(molecule, max_px=max_px, rot_axis='z')

"""
Switch rotation direction
"""

# counterclockwise (ccw)[Default] or clockwise (cw)
# based on rot_axis
# 'x': view from left
# 'y': view from top
# 'z': view into screen
molgif.rot_gif(molecule, max_px=max_px, direction='cw')

"""
Visualize charges
"""

# random charges [-1, 1]
chgs = [-1 + 2 * random.random() for i in molecule]

# manually set the colorbar range (optional)
cb_range = (-1, 1)

# add the charges to atoms object
molecule.set_initial_charges(chgs)

molgif.rot_gif(molecule, max_px=max_px, use_charges=True,
               cb_range=cb_range)
