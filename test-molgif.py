import molgif
import os
import random
import ase.build
import numpy as np

desk = os.path.join(os.path.expanduser('~'), 'Desktop')
os.chdir(desk)

if not os.path.isdir('TESTING_MOLGIF'):
    os.mkdir('TESTING_MOLGIF')
os.chdir('TESTING_MOLGIF')

# load in molecule (ase.Atoms object)
c60 = ase.build.molecule('C60')
c4h4o = ase.build.molecule('C4H4O')
c2h6so = ase.build.molecule('C2H6SO')

# select atoms to animate
atoms = c4h4o
atoms.rotate(60, 'z')

# pixels in largest dimension (width or height)
max_px = 300

"""
Easily create gif with rot_gif function
"""

# create rotating gif with rot_gif function
molgif.rot_gif(atoms, max_px=max_px)

"""
Use smart_rotate to find best viewing angle and add a legend
"""

molgif.rot_gif(atoms, max_px=max_px, smart_rotate=True,
               draw_legend=True)

# can also smart_rotate atoms object using pca in utils
atoms = molgif.utils.smart_rotate_atoms(atoms)

"""
Specify the color of each atom
"""

# can be a string for one color or a list of custom colors
rainbow = ['red', 'orange', 'yellow',
           'green', 'blue', 'violet'] * (len(atoms) // 2)

# list much match number of atoms
rainbow = rainbow[:len(atoms)]

molgif.rot_gif(atoms, max_px=max_px, colors=rainbow)

"""
Use a dictionary to quickly color by atom type
"""

# specify legend to order by atom size
leg_order = 'size'

# default colors will be used for types not specified
molgif.rot_gif(atoms, max_px=max_px, draw_legend=True,
               leg_order=leg_order, colors=dict(C='hotpink'))

"""
Anchor an atom to be at the center of rotation
"""

# define index of atom to anchor
anchor = min(3, len(atoms) - 1)

colors = ['white'] * len(atoms)
colors[anchor] = '#0892d0'

molgif.rot_gif(atoms, max_px=max_px * 1.54, colors=colors,
               anchor=anchor)

"""
Adjust loop time and FPS
"""

# loop_time = time to complete one rotation (seconds)
molgif.rot_gif(atoms, max_px=max_px, loop_time=2, fps=60)

"""
Turn off bonds and scale atomic sizes
"""

molgif.rot_gif(atoms, max_px=max_px, draw_bonds=False, scale=0.9)

"""
Change rotation axis
"""

# switch between x, y (Default), or z
molgif.rot_gif(atoms, max_px=max_px, rot_axis='x')

"""
Switch rotation direction
"""

# "negative" rot_axis = opposite direction
molgif.rot_gif(atoms, max_px=max_px, rot_axis='-x')

"""
Visualize charges
"""

# random charges [-1, 1]
chgs = np.linspace(-1, 1, len(atoms))
np.random.shuffle(chgs)

# manually set the colorbar range (optional)
cb_range = (-1, 1)

# add the charges to atoms object
atoms.set_initial_charges(chgs)

molgif.rot_gif(atoms, max_px=max_px, use_charges=True)
