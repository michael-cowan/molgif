import molgif
import os
import random
import ase.build

desk = os.path.join(os.path.expanduser('~'), 'Desktop')

# load in molecule (ase.Atoms object)
molecule = ase.build.molecule('biphenyl')

# pixels in largest dimension (width or height)
max_px = 300

"""
Use auto_rotate to find the best viewing angle
"""
# specify save path
save_path = os.path.join(desk, 'biphenyl.gif')

# create rotating gif with rot_gif function
molgif.rot_gif(molecule, save_path, max_px=max_px, auto_rotate=True)

"""
Add a legend
"""
save_path = os.path.join(desk, 'biphenyl-legend.gif')

molgif.rot_gif(molecule, save_path, max_px=max_px, add_legend=True)

"""
Specify the color of each atom
"""
save_path = os.path.join(desk, 'biphenyl-rainbow.gif')

# can be a string for one color or a list of custom colors
rainbow = ['red', 'orange', 'yellow',
           'green', 'blue', 'violet'] * 4

# list much match number of atoms
rainbow = rainbow[:len(molecule)]

molgif.rot_gif(molecule, save_path, max_px=max_px, colors=rainbow)

"""
Use a dictionary to quickly color by atom type
"""
save_path = os.path.join(desk, 'biphenyl-hotpink.gif')

# default colors will be used for types not specified
molgif.rot_gif(molecule, save_path, max_px=max_px, colors=dict(C='hotpink'),
               add_legend=True)

"""
Anchor an atom to be at the center of rotation
"""
save_path = os.path.join(desk, 'biphenyl-anchor.gif')

# define index of atom to anchor
anchor = 3

colors = ['white'] * len(molecule)
colors[anchor] = '#0892d0'

molgif.rot_gif(molecule, save_path, max_px=max_px * 1.54, colors=colors,
               anchor=anchor)

"""
Adjust loop time and FPS
"""
save_path = os.path.join(desk, 'biphenyl-2s-looptime.gif')

# loop_time = time to complete one rotation (seconds)
molgif.rot_gif(molecule, save_path, max_px=max_px, loop_time=2, fps=60)

"""
Turn off bonds and scale atomic sizes
"""
save_path = os.path.join(desk, 'biphenyl-no-bonds.gif')

molgif.rot_gif(molecule, save_path, max_px=max_px, add_bonds=False,
               scale=0.9)

"""
Change rotation axis
"""
save_path = os.path.join(desk, 'biphenyl-rotz.gif')

# switch between x, y (Default), or z
molgif.rot_gif(molecule, save_path, max_px=max_px, rot_axis='z')

"""
Switch rotation direction
"""
save_path = os.path.join(desk, 'biphenyl-cw.gif')

# counterclockwise (ccw)[Default] or clockwise (cw)
# based on rot_axis
# 'x': view from left
# 'y': view from top
# 'z': view into screen
molgif.rot_gif(molecule, save_path, max_px=max_px, direction='cw')

"""
Visualize charges
"""
save_path = os.path.join(desk, 'biphenyl-charges.gif')

# random charges [-1, 1]
chgs = [-1 + 2 * random.random() for i in molecule]

# manually set the colorbar range (optional)
cb_range = (-1, 1)

# add the charges to atoms object
molecule.set_initial_charges(chgs)

molgif.rot_gif(molecule, save_path, max_px=max_px, use_charges=True,
               cb_range=cb_range)
