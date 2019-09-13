# molgif

create smooth gifs of rotating molecules

![C60](gifs/C60.gif)

## Installation

    pip install molgif

## Examples

### Easily create gif with rot_gif function

    import molgif
    import ase.build

    # load in molecule (ase.Atoms object)
    c4h4o = ase.build.molecule('c4h4o')
    c4h4o.rotate(60, 'z')

    # create rotating gif
    molgif.rot_gif(atoms)

![c4h4o](gifs/C4H4O.gif)

### Use smart_rotate to find best viewing angle and add a legend

    molgif.rot_gif(c4h4o, smart_rotate=True, add_legend=True)

    # can also smart_rotate atoms object using function in utils
    c4h4o = molgif.utils.smart_rotate_atoms(c4h4o)

![c4h4o-1](gifs/C4H4O-1.gif)

### Specify the color of each atom

    # can be a string for one color or a list of custom colors
    rainbow = ['red', 'orange', 'yellow',
               'green', 'blue', 'violet'] * 2

    # list length much match number of atoms
    rainbow = rainbow[:len(c4h4o)]

    molgif.rot_gif(c4h4o, colors=rainbow)

![c4h4o-2](gifs/C4H4O-2.gif)

### Use a dictionary to quickly color by atom type and add a legend that's ordered by size

    # default colors will be used for types not specified
    # can also order legend by size
    molgif.rot_gif(c4h4o, colors=dict(C='hotpink'),
                   add_legend=True, leg_order='size')

![c4h4o-3](gifs/C4H4O-3.gif)

### Anchor an atom to be at the center of rotation

    # define index of atom to anchor
    anchor = 3

    colors = ['white'] * len(c4h4o)
    colors[anchor] = '#0892d0'

    molgif.rot_gif(c4h4o, colors=colors,
                   anchor=anchor)

![c4h4o-4](gifs/C4H4O-4.gif)

### Adjust loop time and FPS

    # loop_time = time to complete one rotation (seconds)
    molgif.rot_gif(c4h4o, loop_time=2, fps=60)

![c4h4o-5](gifs/C4H4O-5.gif)

### Turn off bonds and scale atomic sizes

    molgif.rot_gif(c4h4o, add_bonds=False,
                   scale=0.9)

![c4h4o-6](gifs/C4H4O-6.gif)

### Change rotation axis

    # switch between x, y (Default), or z
    molgif.rot_gif(c4h4o, rot_axis='x')

![c4h4o-7](gifs/C4H4O-7.gif)

### Switch rotation direction

    # "negative" rot_axis = opposite direction
    molgif.rot_gif(c4h4o, rot_axis='-x')

![c4h4o-8](gifs/C4H4O-8.gif)

### Visualize charges

    # random charges [-1, 1]
    chgs = np.linspace(-1, 1, len(atoms))
    np.random.shuffle(chgs)

    # add the charges to atoms object
    atoms.set_initial_charges(chgs)

    molgif.rot_gif(atoms, max_px=max_px, use_charges=True)

![c4h4o-9](gifs/C4H4O-9.gif)

## Requirements

- ase
- matplotlib
- pillow
- ImageMagick (command line tools must be installed)
