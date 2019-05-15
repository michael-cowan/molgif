# molgif

create smooth gifs of rotating molecules

![C60](gifs/C60.gif)

## Examples

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

    # turn off bonds and scale atomic sizes
    molgif.rot_gif(molecule, save_path, auto_rotate=True,
                   add_bonds=False, scale=0.9)

![biphenyl-no-bonds](gifs/biphenyl-no-bonds.gif)

## Requirements

- ase
- matplotlib
