# molgif

create smooth gifs of rotating molecules

![C60](gifs/C60.gif)

## Examples

### Automatically rotate molecule for better view

    import molgif
    import ase.build

    # load in molecule (ase.Atoms object)
    molecule = ase.build.molecule('biphenyl')

    # specify save path
    save_path = 'biphenyl.gif'

    # create rotating gif with rot_gif function
    # auto_rotate: rotates molecule for better view
    molgif.rot_gif(molecule, save_path, auto_rotate=True)

![biphenyl](gifs/biphenyl.gif)

### Turn off bonds and scale atomic sizes

    molgif.rot_gif(molecule, save_path, auto_rotate=True,
                   add_bonds=False, scale=0.9)

![biphenyl-no-bonds](gifs/biphenyl-no-bonds.gif)

### Visualize charges and include a colorbar

    import random

    # random charges [-1, 1]
    chgs = [-1 + 2 * random.random() for i in molecule]

    molecule.set_initial_charges(chgs)

    molgif.rot_gif(molecule, save_path, auto_rotate=True,
                   use_charges=True)

![biphenyl-charges](gifs/biphenyl-charges.gif)

## Requirements

- ase
- matplotlib
